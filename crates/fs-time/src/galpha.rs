//! Generalized-α (Chung–Hulbert) for structural dynamics
//! M·q̈ + C·q̇ + K·q = f(t): second-order accurate with CONTROLLABLE
//! high-frequency dissipation via ρ∞ ∈ [0, 1] (ρ∞ = 1: no dissipation;
//! ρ∞ = 0: annihilate the highest mode in one step). The spectral
//! behavior is TESTED against theory in the battery, not just cited.

use fs_la::factor::{Lu, lu};
use fs_solver::{
    LinearOp, NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
};
use std::fmt;

/// Prefactored generalized-α stepper for fixed (M, C, K, h).
pub struct GeneralizedAlpha {
    n: usize,
    h: f64,
    alpha_m: f64,
    alpha_f: f64,
    beta: f64,
    gamma: f64,
    m_mat: Vec<f64>,
    c_mat: Vec<f64>,
    k_mat: Vec<f64>,
    eff: Lu,
}

impl GeneralizedAlpha {
    /// Build from mass/damping/stiffness (row-major n×n) and ρ∞.
    ///
    /// # Panics
    /// Structured panic if the effective matrix is singular (a modeling
    /// error: h and the system matrices are incompatible).
    #[must_use]
    pub fn new(
        m_mat: &[f64],
        c_mat: &[f64],
        k_mat: &[f64],
        n: usize,
        h: f64,
        rho_inf: f64,
    ) -> GeneralizedAlpha {
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        // Chung–Hulbert parameterization.
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f) * (1.0 - alpha_m + alpha_f);
        // Effective matrix: (1−αm)/(βh²)·M + (1−αf)γ/(βh)·C + (1−αf)·K.
        let cm = (1.0 - alpha_m) / (beta * h * h);
        let cc = (1.0 - alpha_f) * gamma / (beta * h);
        let ck = 1.0 - alpha_f;
        let mut eff = vec![0.0f64; n * n];
        for i in 0..n * n {
            eff[i] = cm.mul_add(m_mat[i], cc.mul_add(c_mat[i], ck * k_mat[i]));
        }
        let eff = lu(&eff, n).expect("generalized-alpha effective matrix must be nonsingular");
        GeneralizedAlpha {
            n,
            h,
            alpha_m,
            alpha_f,
            beta,
            gamma,
            m_mat: m_mat.to_vec(),
            c_mat: c_mat.to_vec(),
            k_mat: k_mat.to_vec(),
            eff,
        }
    }

    /// The (γ, β) Newmark parameters in use (diagnostics).
    #[must_use]
    pub fn newmark(&self) -> (f64, f64) {
        (self.gamma, self.beta)
    }
}

fn matvec(a: &[f64], n: usize, x: &[f64], out: &mut [f64]) {
    for i in 0..n {
        let mut acc = 0.0f64;
        for j in 0..n {
            acc = a[i * n + j].mul_add(x[j], acc);
        }
        out[i] = acc;
    }
}

/// One generalized-α step: (q, v, a) at t → t+h with load `f_next`
/// evaluated at t + (1−αf)·h by the caller (constant loads just pass
/// the value). Updates in place.
pub fn galpha_step(
    ga: &GeneralizedAlpha,
    q: &mut [f64],
    v: &mut [f64],
    a: &mut [f64],
    f_next: &[f64],
) {
    let n = ga.n;
    let (h, am, af, beta, gamma) = (ga.h, ga.alpha_m, ga.alpha_f, ga.beta, ga.gamma);
    // Predictors (Newmark form).
    let cm = (1.0 - am) / (beta * h * h);
    let cc = (1.0 - af) * gamma / (beta * h);
    // RHS: the equilibrium load `F_{n+1−αf}` (= `f_next`, already the load at the
    // intermediate time t+(1−αf)h per the API) enters with COEFFICIENT 1, plus
    // the M/C history and the −αf·K·q stiffness history:
    // r = f_next + M·[cm·q + (1−am)/(βh)·v + ((1−am)/(2β) − 1)·a]
    //   + C·[cc·q + ((1−αf)γ/β − 1)·v + (1−αf)h·(γ/(2β) − 1)·a]
    //   − αf·K·q.
    let mv_c1 = (1.0 - am) / (beta * h);
    let mv_c2 = (1.0 - am) / (2.0 * beta) - 1.0;
    let cv_c1 = (1.0 - af) * gamma / beta - 1.0;
    let cv_c2 = (1.0 - af) * h * (gamma / (2.0 * beta) - 1.0);
    let mut tm = vec![0.0f64; n];
    let mut tc = vec![0.0f64; n];
    let mut tk = vec![0.0f64; n];
    let mvec: Vec<f64> = (0..n)
        .map(|i| cm.mul_add(q[i], mv_c1.mul_add(v[i], mv_c2 * a[i])))
        .collect();
    let cvec: Vec<f64> = (0..n)
        .map(|i| cc.mul_add(q[i], cv_c1.mul_add(v[i], cv_c2 * a[i])))
        .collect();
    matvec(&ga.m_mat, n, &mvec, &mut tm);
    matvec(&ga.c_mat, n, &cvec, &mut tc);
    matvec(&ga.k_mat, n, q, &mut tk);
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        rhs[i] = f_next[i] + tm[i] + tc[i] - af * tk[i];
    }
    ga.eff.solve(&mut rhs); // rhs now holds q_{n+1}
    // Newmark corrector for a_{n+1}, v_{n+1}.
    for i in 0..n {
        let dq = rhs[i] - q[i];
        let a_new = (dq / (beta * h * h)) - v[i] / (beta * h) - (0.5 / beta - 1.0) * a[i];
        let v_new = v[i] + h * ((1.0 - gamma) * a[i] + gamma * a_new);
        q[i] = rhs[i];
        v[i] = v_new;
        a[i] = a_new;
    }
}

/// Newton--Krylov policy retained by every operator-backed generalized-alpha
/// integrator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitSolveConfig {
    /// Inner and globalization controls supplied to `fs-solver`.
    pub newton: NewtonKrylovConfig,
    /// Maximum outer attempts in one time step.
    pub max_newton_iterations: usize,
}

impl Default for ImplicitSolveConfig {
    fn default() -> Self {
        Self {
            newton: NewtonKrylovConfig::default(),
            max_newton_iterations: 16,
        }
    }
}

/// Complete iteration receipt for one accepted implicit time step.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitStepTelemetry {
    /// Zero-based accepted step index.
    pub step: usize,
    /// Time at the beginning of the step.
    pub t_start: f64,
    /// Step size.
    pub h: f64,
    /// Newton--Krylov report, including per-outer-iteration Krylov counts.
    pub newton: NewtonReport,
}

impl ImplicitStepTelemetry {
    /// Total inner iterations over every Newton attempt in this step.
    #[must_use]
    pub fn krylov_iterations(&self) -> usize {
        self.newton.history.iter().fold(0usize, |total, iteration| {
            total.saturating_add(iteration.linear_iterations)
        })
    }
}

/// Typed refusal from an operator-backed generalized-alpha step.
#[derive(Debug, Clone, PartialEq)]
pub enum TimeSolveError {
    /// An operator, state vector, or forcing vector has the wrong dimension.
    Dimension {
        /// Semantic role of the mismatched object.
        role: &'static str,
        /// Required dimension.
        expected: usize,
        /// Supplied dimension.
        actual: usize,
    },
    /// `fs-solver` refused the initial Newton checkpoint.
    NewtonSetup(NewtonError),
    /// The configured nonlinear budget ended without an accepted solution.
    NotConverged(NewtonReport),
}

impl fmt::Display for TimeSolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dimension {
                role,
                expected,
                actual,
            } => write!(
                f,
                "{role} dimension {actual} differs from integrator dimension {expected}"
            ),
            Self::NewtonSetup(error) => write!(f, "could not initialize Newton step: {error}"),
            Self::NotConverged(report) => write!(
                f,
                "implicit time step did not converge after {} Newton attempts: {:?}",
                report.iterations, report.diagnosis
            ),
        }
    }
}

impl core::error::Error for TimeSolveError {}

/// Structural residual `M a + C v + r(q) = f` used by the generalized-alpha
/// driver. The tangent action must differentiate the exact `internal_force`
/// implementation at the supplied state.
pub trait SecondOrderProblem {
    /// State dimension.
    fn dimension(&self) -> usize;
    /// Overwrite `output` with `M input`.
    fn mass_apply(&self, input: &[f64], output: &mut [f64]);
    /// Overwrite `output` with `C input`.
    fn damping_apply(&self, input: &[f64], output: &mut [f64]);
    /// Overwrite `output` with the internal force `r(q)`.
    fn internal_force(&self, q: &[f64], output: &mut [f64]);
    /// Overwrite `output` with `Dr(q) direction`.
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]);
}

/// Linear `M`, `C`, and `K` adapter over the shared `fs-solver::LinearOp`
/// interface.
pub struct LinearSecondOrderSystem<'a, M: ?Sized, C: ?Sized, K: ?Sized> {
    mass: &'a M,
    damping: &'a C,
    stiffness: &'a K,
    n: usize,
}

impl<'a, M, C, K> LinearSecondOrderSystem<'a, M, C, K>
where
    M: LinearOp + ?Sized,
    C: LinearOp + ?Sized,
    K: LinearOp + ?Sized,
{
    /// Bind three same-sized linear operators.
    #[must_use]
    pub fn new(mass: &'a M, damping: &'a C, stiffness: &'a K) -> Self {
        let n = mass.n();
        assert_eq!(damping.n(), n, "generalized-alpha damping dimension");
        assert_eq!(stiffness.n(), n, "generalized-alpha stiffness dimension");
        Self {
            mass,
            damping,
            stiffness,
            n,
        }
    }
}

impl<M, C, K> SecondOrderProblem for LinearSecondOrderSystem<'_, M, C, K>
where
    M: LinearOp + ?Sized,
    C: LinearOp + ?Sized,
    K: LinearOp + ?Sized,
{
    fn dimension(&self) -> usize {
        self.n
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.mass.apply(input, output);
    }

    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        self.damping.apply(input, output);
    }

    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        self.stiffness.apply(q, output);
    }

    fn tangent_apply(&self, _q: &[f64], direction: &[f64], output: &mut [f64]) {
        self.stiffness.apply(direction, output);
    }
}

/// Plain-data structural trajectory checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct SecondOrderState {
    /// Current time.
    pub t: f64,
    /// Generalized displacement.
    pub q: Vec<f64>,
    /// Generalized velocity.
    pub v: Vec<f64>,
    /// Generalized acceleration.
    pub a: Vec<f64>,
    /// Accepted steps.
    pub steps: usize,
    /// Complete accepted-step telemetry.
    pub history: Vec<ImplicitStepTelemetry>,
}

impl SecondOrderState {
    /// Construct a checkpoint at one exact time.
    #[must_use]
    pub fn new(t: f64, q: &[f64], v: &[f64], a: &[f64]) -> Self {
        assert_eq!(q.len(), v.len(), "generalized-alpha q/v dimension");
        assert_eq!(q.len(), a.len(), "generalized-alpha q/a dimension");
        Self {
            t,
            q: q.to_vec(),
            v: v.to_vec(),
            a: a.to_vec(),
            steps: 0,
            history: Vec::new(),
        }
    }
}

/// Matrix-free/nonlinear Chung--Hulbert generalized-alpha stepper.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OperatorGeneralizedAlpha {
    n: usize,
    h: f64,
    alpha_m: f64,
    alpha_f: f64,
    beta: f64,
    gamma: f64,
    solve: ImplicitSolveConfig,
}

impl OperatorGeneralizedAlpha {
    /// Configure an operator-backed structural integrator.
    #[must_use]
    pub fn new(n: usize, h: f64, rho_inf: f64, solve: ImplicitSolveConfig) -> Self {
        assert!(n > 0, "generalized-alpha dimension must be positive");
        assert!(
            h.is_finite() && h > 0.0,
            "generalized-alpha h must be positive and finite"
        );
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        assert!(
            solve.max_newton_iterations > 0,
            "generalized-alpha Newton budget must be positive"
        );
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f) * (1.0 - alpha_m + alpha_f);
        Self {
            n,
            h,
            alpha_m,
            alpha_f,
            beta,
            gamma,
            solve,
        }
    }

    /// The `(gamma, beta)` Newmark parameters in use.
    #[must_use]
    pub const fn newmark(&self) -> (f64, f64) {
        (self.gamma, self.beta)
    }

    /// Advance one step. `forcing` is evaluated at
    /// `t_n + (1 - alpha_f) h`. State is changed only after convergence.
    pub fn step<P: SecondOrderProblem + ?Sized>(
        &self,
        state: &mut SecondOrderState,
        problem: &P,
        forcing: &[f64],
    ) -> Result<ImplicitStepTelemetry, TimeSolveError> {
        require_dimension("structural problem", self.n, problem.dimension())?;
        require_dimension("displacement state", self.n, state.q.len())?;
        require_dimension("velocity state", self.n, state.v.len())?;
        require_dimension("acceleration state", self.n, state.a.len())?;
        require_dimension("forcing", self.n, forcing.len())?;

        let residual = StructuralStepResidual {
            method: self,
            problem,
            q0: &state.q,
            v0: &state.v,
            a0: &state.a,
            forcing,
        };
        let h2 = self.h * self.h;
        let guess: Vec<f64> = (0..self.n)
            .map(|i| {
                self.h.mul_add(
                    state.v[i],
                    (h2 * (0.5 - self.beta)).mul_add(state.a[i], state.q[i]),
                )
            })
            .collect();
        let mut newton = NewtonKrylovState::new(&residual, guess, self.solve.newton)
            .map_err(TimeSolveError::NewtonSetup)?;
        let report = newton.run(&residual, self.solve.max_newton_iterations);
        if !report.converged {
            return Err(TimeSolveError::NotConverged(report));
        }
        let q_new = newton.x;
        let (v_new, a_new) = structural_kinematics(self, &state.q, &state.v, &state.a, &q_new);
        let telemetry = ImplicitStepTelemetry {
            step: state.steps,
            t_start: state.t,
            h: self.h,
            newton: report,
        };
        state.q = q_new;
        state.v = v_new;
        state.a = a_new;
        state.t += self.h;
        state.steps += 1;
        state.history.push(telemetry.clone());
        Ok(telemetry)
    }
}

struct StructuralStepResidual<'a, P: ?Sized> {
    method: &'a OperatorGeneralizedAlpha,
    problem: &'a P,
    q0: &'a [f64],
    v0: &'a [f64],
    a0: &'a [f64],
    forcing: &'a [f64],
}

impl<P: SecondOrderProblem + ?Sized> NonlinearProblem for StructuralStepResidual<'_, P> {
    fn dimension(&self) -> usize {
        self.method.n
    }

    fn residual(&self, q_new: &[f64], output: &mut [f64]) {
        let method = self.method;
        let (v_new, a_new) = structural_kinematics(method, self.q0, self.v0, self.a0, q_new);
        let mut q_eval = vec![0.0; method.n];
        let mut v_eval = vec![0.0; method.n];
        let mut a_eval = vec![0.0; method.n];
        for i in 0..method.n {
            q_eval[i] = (1.0 - method.alpha_f).mul_add(q_new[i], method.alpha_f * self.q0[i]);
            v_eval[i] = (1.0 - method.alpha_f).mul_add(v_new[i], method.alpha_f * self.v0[i]);
            a_eval[i] = (1.0 - method.alpha_m).mul_add(a_new[i], method.alpha_m * self.a0[i]);
        }
        let mut mass = vec![0.0; method.n];
        let mut damping = vec![0.0; method.n];
        let mut internal = vec![0.0; method.n];
        self.problem.mass_apply(&a_eval, &mut mass);
        self.problem.damping_apply(&v_eval, &mut damping);
        self.problem.internal_force(&q_eval, &mut internal);
        for i in 0..method.n {
            output[i] = mass[i] + damping[i] + internal[i] - self.forcing[i];
        }
    }

    fn jacobian_apply(&self, q_new: &[f64], direction: &[f64], output: &mut [f64]) {
        let method = self.method;
        let mass_scale = (1.0 - method.alpha_m) / (method.beta * method.h * method.h);
        let damping_scale = (1.0 - method.alpha_f) * method.gamma / (method.beta * method.h);
        let mut mass = vec![0.0; method.n];
        let mut damping = vec![0.0; method.n];
        self.problem.mass_apply(direction, &mut mass);
        self.problem.damping_apply(direction, &mut damping);
        let mut q_eval = vec![0.0; method.n];
        let mut tangent_direction = vec![0.0; method.n];
        for i in 0..method.n {
            q_eval[i] = (1.0 - method.alpha_f).mul_add(q_new[i], method.alpha_f * self.q0[i]);
            tangent_direction[i] = (1.0 - method.alpha_f) * direction[i];
        }
        let mut tangent = vec![0.0; method.n];
        self.problem
            .tangent_apply(&q_eval, &tangent_direction, &mut tangent);
        for i in 0..method.n {
            output[i] = mass_scale.mul_add(mass[i], damping_scale.mul_add(damping[i], tangent[i]));
        }
    }
}

fn structural_kinematics(
    method: &OperatorGeneralizedAlpha,
    q0: &[f64],
    v0: &[f64],
    a0: &[f64],
    q_new: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let mut v_new = vec![0.0; method.n];
    let mut a_new = vec![0.0; method.n];
    for i in 0..method.n {
        let dq = q_new[i] - q0[i];
        a_new[i] = dq / (method.beta * method.h * method.h)
            - v0[i] / (method.beta * method.h)
            - (0.5 / method.beta - 1.0) * a0[i];
        v_new[i] = v0[i] + method.h * ((1.0 - method.gamma) * a0[i] + method.gamma * a_new[i]);
    }
    (v_new, a_new)
}

/// Prefactored first-order generalized-alpha fast lane for
/// `M udot + A u = f`.
pub struct FirstOrderGeneralizedAlpha {
    n: usize,
    h: f64,
    alpha_m: f64,
    alpha_f: f64,
    gamma: f64,
    mass: Vec<f64>,
    operator: Vec<f64>,
    effective: Lu,
}

impl FirstOrderGeneralizedAlpha {
    /// Build a dense first-order stepper from row-major `M` and `A`.
    #[must_use]
    pub fn new(mass: &[f64], operator: &[f64], n: usize, h: f64, rho_inf: f64) -> Self {
        assert!(
            n > 0,
            "first-order generalized-alpha dimension must be positive"
        );
        assert_eq!(
            mass.len(),
            n * n,
            "first-order generalized-alpha mass shape"
        );
        assert_eq!(
            operator.len(),
            n * n,
            "first-order generalized-alpha operator shape"
        );
        assert!(h.is_finite() && h > 0.0, "first-order generalized-alpha h");
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        let (alpha_m, alpha_f, gamma) = first_order_parameters(rho_inf);
        let mass_scale = alpha_m / (gamma * h);
        let mut effective = vec![0.0; n * n];
        for i in 0..n * n {
            effective[i] = mass_scale.mul_add(mass[i], alpha_f * operator[i]);
        }
        let effective = lu(&effective, n)
            .expect("first-order generalized-alpha effective matrix must be nonsingular");
        Self {
            n,
            h,
            alpha_m,
            alpha_f,
            gamma,
            mass: mass.to_vec(),
            operator: operator.to_vec(),
            effective,
        }
    }
}

/// One dense first-order generalized-alpha step. `forcing` is evaluated at
/// `t_n + alpha_f h` by the caller.
pub fn first_order_galpha_step(
    method: &FirstOrderGeneralizedAlpha,
    u: &mut [f64],
    rate: &mut [f64],
    forcing: &[f64],
) {
    let n = method.n;
    assert_eq!(u.len(), n, "first-order generalized-alpha state dimension");
    assert_eq!(
        rate.len(),
        n,
        "first-order generalized-alpha rate dimension"
    );
    assert_eq!(
        forcing.len(),
        n,
        "first-order generalized-alpha forcing dimension"
    );
    let mass_scale = method.alpha_m / (method.gamma * method.h);
    let mut rate_history = vec![0.0; n];
    let mut state_history = vec![0.0; n];
    for i in 0..n {
        rate_history[i] =
            mass_scale.mul_add(u[i], -(1.0 - method.alpha_m / method.gamma) * rate[i]);
        state_history[i] = (1.0 - method.alpha_f) * u[i];
    }
    let mut mass_history = vec![0.0; n];
    let mut operator_history = vec![0.0; n];
    matvec(&method.mass, n, &rate_history, &mut mass_history);
    matvec(&method.operator, n, &state_history, &mut operator_history);
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        rhs[i] = forcing[i] + mass_history[i] - operator_history[i];
    }
    method.effective.solve(&mut rhs);
    for i in 0..n {
        rate[i] = (rhs[i] - u[i]) / (method.gamma * method.h)
            - ((1.0 - method.gamma) / method.gamma) * rate[i];
        u[i] = rhs[i];
    }
}

/// First-order residual `M udot + r(t, u) = f`.
pub trait FirstOrderProblem {
    /// State dimension.
    fn dimension(&self) -> usize;
    /// Overwrite `output` with `M input`.
    fn mass_apply(&self, input: &[f64], output: &mut [f64]);
    /// Overwrite `output` with `r(t, u)`.
    fn internal_force(&self, t: f64, u: &[f64], output: &mut [f64]);
    /// Overwrite `output` with `Dr(t, u) direction`.
    fn tangent_apply(&self, t: f64, u: &[f64], direction: &[f64], output: &mut [f64]);
}

/// Linear first-order adapter over `fs-solver::LinearOp`.
pub struct LinearFirstOrderSystem<'a, M: ?Sized, A: ?Sized> {
    mass: &'a M,
    operator: &'a A,
    n: usize,
}

impl<'a, M, A> LinearFirstOrderSystem<'a, M, A>
where
    M: LinearOp + ?Sized,
    A: LinearOp + ?Sized,
{
    /// Bind same-sized mass and evolution operators.
    #[must_use]
    pub fn new(mass: &'a M, operator: &'a A) -> Self {
        let n = mass.n();
        assert_eq!(operator.n(), n, "first-order generalized-alpha dimension");
        Self { mass, operator, n }
    }
}

impl<M, A> FirstOrderProblem for LinearFirstOrderSystem<'_, M, A>
where
    M: LinearOp + ?Sized,
    A: LinearOp + ?Sized,
{
    fn dimension(&self) -> usize {
        self.n
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.mass.apply(input, output);
    }

    fn internal_force(&self, _t: f64, u: &[f64], output: &mut [f64]) {
        self.operator.apply(u, output);
    }

    fn tangent_apply(&self, _t: f64, _u: &[f64], direction: &[f64], output: &mut [f64]) {
        self.operator.apply(direction, output);
    }
}

/// Plain-data first-order generalized-alpha trajectory checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct FirstOrderState {
    /// Current time.
    pub t: f64,
    /// Current state.
    pub u: Vec<f64>,
    /// Current time derivative.
    pub rate: Vec<f64>,
    /// Accepted steps.
    pub steps: usize,
    /// Complete accepted-step telemetry.
    pub history: Vec<ImplicitStepTelemetry>,
}

impl FirstOrderState {
    /// Construct a checkpoint from a state and a consistent initial rate.
    #[must_use]
    pub fn new(t: f64, u: &[f64], rate: &[f64]) -> Self {
        assert_eq!(
            u.len(),
            rate.len(),
            "first-order generalized-alpha state/rate dimension"
        );
        Self {
            t,
            u: u.to_vec(),
            rate: rate.to_vec(),
            steps: 0,
            history: Vec::new(),
        }
    }
}

/// Operator-backed first-order generalized-alpha method (Jansen--Whiting--
/// Hulbert parameterization).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OperatorFirstOrderGeneralizedAlpha {
    n: usize,
    h: f64,
    alpha_m: f64,
    alpha_f: f64,
    gamma: f64,
    solve: ImplicitSolveConfig,
}

impl OperatorFirstOrderGeneralizedAlpha {
    /// Configure a first-order operator/nonlinear stepper.
    #[must_use]
    pub fn new(n: usize, h: f64, rho_inf: f64, solve: ImplicitSolveConfig) -> Self {
        assert!(
            n > 0,
            "first-order generalized-alpha dimension must be positive"
        );
        assert!(h.is_finite() && h > 0.0, "first-order generalized-alpha h");
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        assert!(
            solve.max_newton_iterations > 0,
            "first-order generalized-alpha Newton budget must be positive"
        );
        let (alpha_m, alpha_f, gamma) = first_order_parameters(rho_inf);
        Self {
            n,
            h,
            alpha_m,
            alpha_f,
            gamma,
            solve,
        }
    }

    /// Advance one step. `forcing` is evaluated at `t_n + alpha_f h`.
    /// State is changed only after convergence.
    pub fn step<P: FirstOrderProblem + ?Sized>(
        &self,
        state: &mut FirstOrderState,
        problem: &P,
        forcing: &[f64],
    ) -> Result<ImplicitStepTelemetry, TimeSolveError> {
        require_dimension("first-order problem", self.n, problem.dimension())?;
        require_dimension("first-order state", self.n, state.u.len())?;
        require_dimension("first-order rate", self.n, state.rate.len())?;
        require_dimension("forcing", self.n, forcing.len())?;
        let residual = FirstOrderStepResidual {
            method: self,
            problem,
            t0: state.t,
            u0: &state.u,
            rate0: &state.rate,
            forcing,
        };
        let guess: Vec<f64> = state
            .u
            .iter()
            .zip(&state.rate)
            .map(|(u, rate)| self.h.mul_add(*rate, *u))
            .collect();
        let mut newton = NewtonKrylovState::new(&residual, guess, self.solve.newton)
            .map_err(TimeSolveError::NewtonSetup)?;
        let report = newton.run(&residual, self.solve.max_newton_iterations);
        if !report.converged {
            return Err(TimeSolveError::NotConverged(report));
        }
        let u_new = newton.x;
        let rate_new = first_order_rate(self, &state.u, &state.rate, &u_new);
        let telemetry = ImplicitStepTelemetry {
            step: state.steps,
            t_start: state.t,
            h: self.h,
            newton: report,
        };
        state.u = u_new;
        state.rate = rate_new;
        state.t += self.h;
        state.steps += 1;
        state.history.push(telemetry.clone());
        Ok(telemetry)
    }
}

struct FirstOrderStepResidual<'a, P: ?Sized> {
    method: &'a OperatorFirstOrderGeneralizedAlpha,
    problem: &'a P,
    t0: f64,
    u0: &'a [f64],
    rate0: &'a [f64],
    forcing: &'a [f64],
}

impl<P: FirstOrderProblem + ?Sized> NonlinearProblem for FirstOrderStepResidual<'_, P> {
    fn dimension(&self) -> usize {
        self.method.n
    }

    fn residual(&self, u_new: &[f64], output: &mut [f64]) {
        let method = self.method;
        let rate_new = first_order_rate(method, self.u0, self.rate0, u_new);
        let mut u_eval = vec![0.0; method.n];
        let mut rate_eval = vec![0.0; method.n];
        for i in 0..method.n {
            u_eval[i] = method
                .alpha_f
                .mul_add(u_new[i], (1.0 - method.alpha_f) * self.u0[i]);
            rate_eval[i] = method
                .alpha_m
                .mul_add(rate_new[i], (1.0 - method.alpha_m) * self.rate0[i]);
        }
        let t_eval = method.h.mul_add(method.alpha_f, self.t0);
        let mut mass = vec![0.0; method.n];
        let mut internal = vec![0.0; method.n];
        self.problem.mass_apply(&rate_eval, &mut mass);
        self.problem.internal_force(t_eval, &u_eval, &mut internal);
        for i in 0..method.n {
            output[i] = mass[i] + internal[i] - self.forcing[i];
        }
    }

    fn jacobian_apply(&self, u_new: &[f64], direction: &[f64], output: &mut [f64]) {
        let method = self.method;
        let mass_scale = method.alpha_m / (method.gamma * method.h);
        let mut mass = vec![0.0; method.n];
        self.problem.mass_apply(direction, &mut mass);
        let mut u_eval = vec![0.0; method.n];
        let mut tangent_direction = vec![0.0; method.n];
        for i in 0..method.n {
            u_eval[i] = method
                .alpha_f
                .mul_add(u_new[i], (1.0 - method.alpha_f) * self.u0[i]);
            tangent_direction[i] = method.alpha_f * direction[i];
        }
        let t_eval = method.h.mul_add(method.alpha_f, self.t0);
        let mut tangent = vec![0.0; method.n];
        self.problem
            .tangent_apply(t_eval, &u_eval, &tangent_direction, &mut tangent);
        for i in 0..method.n {
            output[i] = mass_scale.mul_add(mass[i], tangent[i]);
        }
    }
}

fn first_order_parameters(rho_inf: f64) -> (f64, f64, f64) {
    let alpha_m = (3.0 - rho_inf) / (2.0 * (1.0 + rho_inf));
    let alpha_f = 1.0 / (1.0 + rho_inf);
    let gamma = 0.5 + alpha_m - alpha_f;
    (alpha_m, alpha_f, gamma)
}

fn first_order_rate(
    method: &OperatorFirstOrderGeneralizedAlpha,
    u0: &[f64],
    rate0: &[f64],
    u_new: &[f64],
) -> Vec<f64> {
    (0..method.n)
        .map(|i| {
            (u_new[i] - u0[i]) / (method.gamma * method.h)
                - ((1.0 - method.gamma) / method.gamma) * rate0[i]
        })
        .collect()
}

fn require_dimension(
    role: &'static str,
    expected: usize,
    actual: usize,
) -> Result<(), TimeSolveError> {
    if expected == actual {
        Ok(())
    } else {
        Err(TimeSolveError::Dimension {
            role,
            expected,
            actual,
        })
    }
}
