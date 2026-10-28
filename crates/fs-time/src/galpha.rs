//! Generalized-α (Chung–Hulbert) for structural dynamics
//! M·q̈ + C·q̇ + K·q = f(t): second-order accurate with CONTROLLABLE
//! high-frequency dissipation via ρ∞ ∈ [0, 1] (ρ∞ = 1: no dissipation;
//! ρ∞ = 0: annihilate the highest mode in one step). The spectral
//! behavior is TESTED against theory in the battery, not just cited.

use fs_la::factor::{Lu, lu};
use fs_solver::{
    CancellableSolveReport, LinearOp, NewtonError, NewtonKrylovConfig, NewtonKrylovState,
    NewtonReport, NonlinearProblem, SolverCallback, SolverRunProgress,
};
use std::fmt;

#[path = "galpha/adjoint.rs"]
pub mod adjoint;
#[path = "galpha/initialization.rs"]
pub mod initialization;
#[path = "galpha/second_order_adjoint.rs"]
pub mod second_order_adjoint;

/// Prefactored generalized-α stepper for fixed (M, C, K, h).
pub struct GeneralizedAlpha {
    n: usize,
    h: f64,
    beta: f64,
    gamma: f64,
    k_mat: Vec<f64>,
    velocity_history: Vec<f64>,
    acceleration_history: Vec<f64>,
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
        assert!(n > 0, "generalized-alpha dimension must be positive");
        let matrix_len = n
            .checked_mul(n)
            .expect("generalized-alpha matrix dimension overflows usize");
        assert_eq!(m_mat.len(), matrix_len, "generalized-alpha mass shape");
        assert_eq!(c_mat.len(), matrix_len, "generalized-alpha damping shape");
        assert_eq!(k_mat.len(), matrix_len, "generalized-alpha stiffness shape");
        assert!(
            h.is_finite() && h > 0.0,
            "generalized-alpha h must be positive and finite"
        );
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        // Chung–Hulbert parameterization.
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f) * (1.0 - alpha_m + alpha_f);
        // Solve directly for endpoint acceleration. The displacement form
        // magnifies roundoff by 1/h² when recovering acceleration from two
        // absolute positions; its M*q/h² history also overwhelms small forces.
        // A_eff = (1−αm)M + (1−αf)γh C + (1−αf)βh² K.
        let cm = 1.0 - alpha_m;
        let cc = (1.0 - alpha_f) * gamma * h;
        let ck = (1.0 - alpha_f) * beta * h * h;
        let kv = (1.0 - alpha_f) * h;
        let ca = (1.0 - alpha_f) * (1.0 - gamma) * h;
        let ka = (1.0 - alpha_f) * (0.5 - beta) * h * h;
        let mut eff = vec![0.0f64; matrix_len];
        let mut velocity_history = vec![0.0; matrix_len];
        let mut acceleration_history = vec![0.0; matrix_len];
        for i in 0..matrix_len {
            eff[i] = cm.mul_add(m_mat[i], cc.mul_add(c_mat[i], ck * k_mat[i]));
            velocity_history[i] = kv.mul_add(k_mat[i], c_mat[i]);
            acceleration_history[i] =
                alpha_m.mul_add(m_mat[i], ca.mul_add(c_mat[i], ka * k_mat[i]));
        }
        let eff = lu(&eff, n).expect("generalized-alpha effective matrix must be nonsingular");
        GeneralizedAlpha {
            n,
            h,
            beta,
            gamma,
            k_mat: k_mat.to_vec(),
            velocity_history,
            acceleration_history,
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
    assert_eq!(q.len(), n, "generalized-alpha displacement state dimension");
    assert_eq!(v.len(), n, "generalized-alpha velocity state dimension");
    assert_eq!(a.len(), n, "generalized-alpha acceleration state dimension");
    assert_eq!(f_next.len(), n, "generalized-alpha forcing dimension");
    let (h, beta, gamma) = (ga.h, ga.beta, ga.gamma);
    // Eliminating the Newmark predictors gives
    // r = f_next - K*q - [C + (1−αf)h K]*v
    //     - [αm M + (1−αf)(1−γ)h C + (1−αf)(½−β)h² K]*a.
    // The two history operators are assembled once with the factorization.
    // History assembly uses one n-vector instead of six; the unchanged LU
    // solve additionally allocates its permutation scratch.
    let mut rhs = f_next.to_vec();
    for i in 0..n {
        for j in 0..n {
            let entry = i * n + j;
            rhs[i] = (-ga.k_mat[entry]).mul_add(q[j], rhs[i]);
            rhs[i] = (-ga.velocity_history[entry]).mul_add(v[j], rhs[i]);
            rhs[i] = (-ga.acceleration_history[entry]).mul_add(a[j], rhs[i]);
        }
    }
    ga.eff.solve(&mut rhs); // rhs now holds a_{n+1}
    // Newmark increments retain velocity and acceleration even when the
    // position increment rounds away relative to a prestressed equilibrium.
    for i in 0..n {
        let dq = h.mul_add(v[i], (h * h) * beta.mul_add(rhs[i], (0.5 - beta) * a[i]));
        let dv = h * gamma.mul_add(rhs[i], (1.0 - gamma) * a[i]);
        q[i] += dq;
        v[i] += dv;
        a[i] = rhs[i];
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
    /// Cooperative cancellation occurred before the public state was changed.
    Cancelled,
    /// An operator, state vector, or forcing vector has the wrong dimension.
    Dimension {
        /// Semantic role of the mismatched object.
        role: &'static str,
        /// Required dimension.
        expected: usize,
        /// Supplied dimension.
        actual: usize,
    },
    /// The configured positive step cannot advance the checkpoint clock to a
    /// distinct finite endpoint.
    InvalidTimeAdvance {
        /// Exact source-time bits.
        t_bits: u64,
        /// Exact configured-step bits.
        h_bits: u64,
        /// Exact rounded endpoint bits.
        next_t_bits: u64,
    },
    /// A time-indexed internal stage cannot be represented strictly after the
    /// source clock even though the endpoint itself is valid.
    InvalidStageTime {
        /// Exact source-time bits.
        t_bits: u64,
        /// Exact configured-step bits.
        h_bits: u64,
        /// Exact stage-fraction bits.
        fraction_bits: u64,
        /// Exact rounded stage-time bits.
        stage_t_bits: u64,
    },
    /// The accepted-step counter cannot represent another commit.
    StepCounterOverflow {
        /// Exact current counter value.
        steps: usize,
    },
    /// The physical model refused the actual step or a state before publication.
    ModelRefused {
        /// Admission point: timestep, initial state, implicit stage or endpoint.
        phase: &'static str,
        /// Original model-owner diagnostic.
        reason: String,
    },
    /// `fs-solver` refused the initial Newton checkpoint.
    NewtonSetup(NewtonError),
    /// The configured nonlinear budget ended without an accepted solution.
    NotConverged(NewtonReport),
    /// A numerical callback unwound during Newton; the entire step was discarded.
    SolverCallbackPanicked(SolverCallback),
}

impl fmt::Display for TimeSolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cancelled => write!(f, "generalized-alpha step cancelled before publication"),
            Self::Dimension {
                role,
                expected,
                actual,
            } => write!(
                f,
                "{role} dimension {actual} differs from integrator dimension {expected}"
            ),
            Self::InvalidTimeAdvance {
                t_bits,
                h_bits,
                next_t_bits,
            } => write!(
                f,
                "generalized-alpha clock cannot advance finitely (t bits 0x{t_bits:016x}, h bits \
                 0x{h_bits:016x}, next bits 0x{next_t_bits:016x})"
            ),
            Self::InvalidStageTime {
                t_bits,
                h_bits,
                fraction_bits,
                stage_t_bits,
            } => write!(
                f,
                "generalized-alpha stage time cannot advance finitely (t bits 0x{t_bits:016x}, \
                 h bits 0x{h_bits:016x}, fraction bits 0x{fraction_bits:016x}, stage bits \
                 0x{stage_t_bits:016x})"
            ),
            Self::StepCounterOverflow { steps } => write!(
                f,
                "generalized-alpha accepted-step counter cannot advance past {steps}"
            ),
            Self::ModelRefused { phase, reason } => {
                write!(f, "structural model refused {phase}: {reason}")
            }
            Self::NewtonSetup(error) => write!(f, "could not initialize Newton step: {error}"),
            Self::NotConverged(report) => write!(
                f,
                "implicit time step did not converge after {} Newton attempts: {:?}",
                report.iterations, report.diagnosis
            ),
            Self::SolverCallbackPanicked(callback) => write!(
                f,
                "generalized-alpha Newton {callback:?} callback panicked; step discarded"
            ),
        }
    }
}

impl core::error::Error for TimeSolveError {}

fn completed_newton(
    outcome: CancellableSolveReport<NewtonReport>,
) -> Result<NewtonReport, TimeSolveError> {
    match outcome.progress {
        SolverRunProgress::Complete => Ok(outcome.report),
        SolverRunProgress::Paused => Err(TimeSolveError::Cancelled),
        SolverRunProgress::CallbackPanicked(callback) => {
            Err(TimeSolveError::SolverCallbackPanicked(callback))
        }
        SolverRunProgress::DimensionMismatch { expected, actual } => {
            Err(TimeSolveError::Dimension {
                role: "Newton problem",
                expected,
                actual,
            })
        }
    }
}

fn checked_time_advance(t: f64, h: f64) -> Result<f64, TimeSolveError> {
    let next_t = t + h;
    if t.is_finite() && next_t.is_finite() && next_t > t {
        Ok(next_t)
    } else {
        Err(TimeSolveError::InvalidTimeAdvance {
            t_bits: t.to_bits(),
            h_bits: h.to_bits(),
            next_t_bits: next_t.to_bits(),
        })
    }
}

/// Coefficients of the structural step's effective Jacobian
/// `mass*M + damping*R_v(q,v) + tangent*R_q(q,v)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondOrderOperatorWeights {
    /// `(1 - alpha_m) / (beta*h*h)`.
    pub mass: f64,
    /// `(1 - alpha_f)*gamma / (beta*h)`.
    pub damping: f64,
    /// `1 - alpha_f`.
    pub tangent: f64,
}

/// Structural residual `M a + R(q, v) = f` used by generalized-alpha.
/// Mass is state independent. By default `R(q,v) = C v + r(q)`, preserving
/// existing implementations. State-dependent resistance must override both
/// partial actions together with `resistance`; all callbacks are pure and
/// differentiate the same force evaluated at the supplied stage state.
pub trait SecondOrderProblem {
    /// Admit the actual integrator timestep in seconds before any solve.
    /// The default imposes no additional model-specific operating limit.
    /// Implementations must be pure, bounded and reproducible during replay.
    fn validate_step(&self, _step_s: f64) -> Result<(), String> {
        Ok(())
    }

    /// Admit physical displacement and velocity before evaluating a stage,
    /// initializing acceleration, or publishing an accepted endpoint.
    /// This restricts the physical domain; it is not a projection or clipping
    /// operation. Adjoint derivatives apply inside that admitted domain.
    fn validate_state(&self, _q: &[f64], _v: &[f64]) -> Result<(), String> {
        Ok(())
    }

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

    /// Overwrite with the complete resisting force R(q,v).
    /// The default combines the legacy linear damping and internal force.
    fn resistance(&self, q: &[f64], v: &[f64], output: &mut [f64]) {
        let mut internal = vec![f64::NAN; self.dimension()];
        self.internal_force(q, &mut internal);
        self.damping_apply(v, output);
        for (value, force) in output.iter_mut().zip(internal) {
            *value += force;
        }
    }

    /// Overwrite with the partial R_q(q,v) direction, holding velocity fixed.
    fn resistance_displacement_apply(
        &self,
        q: &[f64],
        _v: &[f64],
        direction: &[f64],
        output: &mut [f64],
    ) {
        self.tangent_apply(q, direction, output);
    }

    /// Overwrite with the partial R_v(q,v) direction, holding displacement fixed.
    fn resistance_velocity_apply(
        &self,
        _q: &[f64],
        _v: &[f64],
        direction: &[f64],
        output: &mut [f64],
    ) {
        self.damping_apply(direction, output);
    }

    /// Approximately solve the effective Jacobian system for the residual.
    ///
    /// `q` is the same interpolated stage displacement used by `tangent_apply`,
    /// and `weights` give the exact effective mass, damping, and tangent
    /// coefficients. This is a right preconditioner inside flexible GMRES;
    /// `outer_iteration` and `inner_iteration` are logical Newton and Krylov
    /// indices. Overwrite every output component with a finite value. Invalid
    /// output causes a solver refusal before the time-step state is published.
    ///
    /// The default copies the residual (identity). Implementations must bound
    /// their own work and reproduce their results during checkpoint replay.
    /// This hook changes only the forward implicit solve; adjoints continue
    /// to use their separately supplied transpose-solve preconditioner.
    /// It is an approximation: nonlinear resistance may add velocity-dependent
    /// terms beyond the legacy mass/damping/tangent weights.
    fn preconditioner_apply(
        &self,
        _q: &[f64],
        _weights: SecondOrderOperatorWeights,
        _outer_iteration: usize,
        _inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.copy_from_slice(residual);
    }
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
        self.step_controlled(state, problem, forcing, &mut || false)
    }

    /// Advance with cancellation checks before setup, before/after Newton
    /// residual, Jacobian and preconditioner callbacks, within inner Krylov
    /// and backtracking work, and before publication. Individual callbacks
    /// must bound their own work.
    /// Cancellation leaves displacement, velocity, acceleration and history
    /// unchanged, including when a previous step has already been accepted.
    pub fn step_controlled<P, Cancel>(
        &self,
        state: &mut SecondOrderState,
        problem: &P,
        forcing: &[f64],
        cancelled: &mut Cancel,
    ) -> Result<ImplicitStepTelemetry, TimeSolveError>
    where
        P: SecondOrderProblem + ?Sized,
        Cancel: FnMut() -> bool,
    {
        structural_poll(cancelled)?;
        require_dimension("structural problem", self.n, problem.dimension())?;
        require_dimension("displacement state", self.n, state.q.len())?;
        require_dimension("velocity state", self.n, state.v.len())?;
        require_dimension("acceleration state", self.n, state.a.len())?;
        require_dimension("forcing", self.n, forcing.len())?;
        let next_t = checked_time_advance(state.t, self.h)?;
        let next_steps = state
            .steps
            .checked_add(1)
            .ok_or(TimeSolveError::StepCounterOverflow { steps: state.steps })?;

        let admission = problem.validate_step(self.h);
        structural_poll(cancelled)?;
        structural_admit(admission, "timestep")?;
        let admission = problem.validate_state(&state.q, &state.v);
        structural_poll(cancelled)?;
        structural_admit(admission, "initial state")?;

        let residual = StructuralStepResidual {
            method: self,
            problem,
            q0: &state.q,
            v0: &state.v,
            a0: &state.a,
            forcing,
            admission_error: std::cell::RefCell::new(None),
        };
        // Solve the Newmark displacement correction beta*h^2*a_next,
        // rather than absolute displacement. Recovering acceleration from
        // q_next - q_old destroys resolved force/velocity increments when the
        // displacement change is small relative to a prestressed equilibrium.
        // The effective Jacobian and its preconditioner weights are unchanged.
        let guess = vec![0.0; self.n];
        let setup = NewtonKrylovState::new(&residual, guess, self.solve.newton);
        structural_poll(cancelled)?;
        let mut newton = setup.map_err(|error| {
            residual
                .take_admission_error()
                .unwrap_or(TimeSolveError::NewtonSetup(error))
        })?;
        let mut report = newton.run(&residual, 0);
        for _ in 0..self.solve.max_newton_iterations {
            // Plateau is diagnostic, not terminal in the shared solver.
            if report.converged
                || !matches!(
                    report.diagnosis,
                    Some(
                        fs_solver::NewtonStallDiagnosis::BudgetExhausted
                            | fs_solver::NewtonStallDiagnosis::Plateau
                    )
                )
            {
                break;
            }
            structural_poll(cancelled)?;
            report = completed_newton(newton.run_controlled(&residual, 1, cancelled))?;
            structural_poll(cancelled)?;
        }
        if !report.converged {
            return Err(residual
                .take_admission_error()
                .unwrap_or(TimeSolveError::NotConverged(report)));
        }
        let mut q_new = newton.x;
        let mut v_new = vec![0.0; self.n];
        let mut a_new = vec![0.0; self.n];
        for i in 0..self.n {
            let (dq, dv, acceleration) =
                structural_increment(self, state.v[i], state.a[i], q_new[i]);
            q_new[i] = state.q[i] + dq;
            v_new[i] = state.v[i] + dv;
            a_new[i] = acceleration;
        }
        let admission = problem.validate_state(&q_new, &v_new);
        structural_poll(cancelled)?;
        structural_admit(admission, "accepted endpoint")?;
        let telemetry = ImplicitStepTelemetry {
            step: state.steps,
            t_start: state.t,
            h: self.h,
            newton: report,
        };
        structural_poll(cancelled)?;
        state.q = q_new;
        state.v = v_new;
        state.a = a_new;
        state.t = next_t;
        state.steps = next_steps;
        state.history.push(telemetry.clone());
        Ok(telemetry)
    }
}

fn structural_admit(result: Result<(), String>, phase: &'static str) -> Result<(), TimeSolveError> {
    result.map_err(|reason| TimeSolveError::ModelRefused { phase, reason })
}

fn structural_poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), TimeSolveError> {
    if cancelled() {
        Err(TimeSolveError::Cancelled)
    } else {
        Ok(())
    }
}

struct StructuralStepResidual<'a, P: ?Sized> {
    method: &'a OperatorGeneralizedAlpha,
    problem: &'a P,
    q0: &'a [f64],
    v0: &'a [f64],
    a0: &'a [f64],
    forcing: &'a [f64],
    admission_error: std::cell::RefCell<Option<String>>,
}

impl<P: SecondOrderProblem + ?Sized> StructuralStepResidual<'_, P> {
    fn admit_stage(&self, q: &[f64], v: &[f64], output: &mut [f64]) -> bool {
        match self.problem.validate_state(q, v) {
            // A rejected line-search trial may recover. Retain its diagnostic
            // only for a later failed solve; convergence ignores this scratch.
            Ok(()) => true,
            Err(reason) => {
                self.admission_error.replace(Some(reason));
                output.fill(f64::NAN);
                false
            }
        }
    }

    fn take_admission_error(&self) -> Option<TimeSolveError> {
        self.admission_error
            .borrow_mut()
            .take()
            .map(|reason| TimeSolveError::ModelRefused {
                phase: "implicit stage",
                reason,
            })
    }
}

impl<P: SecondOrderProblem + ?Sized> NonlinearProblem for StructuralStepResidual<'_, P> {
    fn dimension(&self) -> usize {
        self.method.n
    }

    fn residual(&self, correction: &[f64], output: &mut [f64]) {
        let method = self.method;
        let mut q_eval = vec![0.0; method.n];
        let mut v_eval = vec![0.0; method.n];
        let mut a_eval = vec![0.0; method.n];
        for i in 0..method.n {
            let (dq, dv, acceleration) =
                structural_increment(method, self.v0[i], self.a0[i], correction[i]);
            q_eval[i] = (1.0 - method.alpha_f).mul_add(dq, self.q0[i]);
            v_eval[i] = (1.0 - method.alpha_f).mul_add(dv, self.v0[i]);
            a_eval[i] = (1.0 - method.alpha_m).mul_add(acceleration, method.alpha_m * self.a0[i]);
        }
        let mut mass = vec![0.0; method.n];
        let mut resistance = vec![0.0; method.n];
        if !self.admit_stage(&q_eval, &v_eval, output) {
            return;
        }
        self.problem.mass_apply(&a_eval, &mut mass);
        self.problem.resistance(&q_eval, &v_eval, &mut resistance);
        for i in 0..method.n {
            output[i] = mass[i] + resistance[i] - self.forcing[i];
        }
    }

    fn jacobian_apply(&self, correction: &[f64], direction: &[f64], output: &mut [f64]) {
        let method = self.method;
        let mass_scale = (1.0 - method.alpha_m) / (method.beta * method.h * method.h);
        let damping_scale = (1.0 - method.alpha_f) * method.gamma / (method.beta * method.h);
        let mut mass = vec![0.0; method.n];
        let mut damping = vec![0.0; method.n];
        self.problem.mass_apply(direction, &mut mass);
        let mut q_eval = vec![0.0; method.n];
        let mut v_eval = vec![0.0; method.n];
        let mut tangent_direction = vec![0.0; method.n];
        for i in 0..method.n {
            let (dq, dv, _) = structural_increment(method, self.v0[i], self.a0[i], correction[i]);
            q_eval[i] = (1.0 - method.alpha_f).mul_add(dq, self.q0[i]);
            v_eval[i] = (1.0 - method.alpha_f).mul_add(dv, self.v0[i]);
            tangent_direction[i] = (1.0 - method.alpha_f) * direction[i];
        }
        if !self.admit_stage(&q_eval, &v_eval, output) {
            return;
        }
        self.problem
            .resistance_velocity_apply(&q_eval, &v_eval, direction, &mut damping);
        let mut tangent = vec![0.0; method.n];
        self.problem.resistance_displacement_apply(
            &q_eval,
            &v_eval,
            &tangent_direction,
            &mut tangent,
        );
        for i in 0..method.n {
            output[i] = mass_scale.mul_add(mass[i], damping_scale.mul_add(damping[i], tangent[i]));
        }
    }

    fn preconditioner_apply(
        &self,
        correction: &[f64],
        outer_iteration: usize,
        inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        let method = self.method;
        let q_eval: Vec<f64> = correction
            .iter()
            .enumerate()
            .map(|(i, &correction)| {
                let (dq, _, _) = structural_increment(method, self.v0[i], self.a0[i], correction);
                (1.0 - method.alpha_f).mul_add(dq, self.q0[i])
            })
            .collect();
        let weights = SecondOrderOperatorWeights {
            mass: (1.0 - method.alpha_m) / (method.beta * method.h * method.h),
            damping: (1.0 - method.alpha_f) * method.gamma / (method.beta * method.h),
            tangent: 1.0 - method.alpha_f,
        };
        self.problem.preconditioner_apply(
            &q_eval,
            weights,
            outer_iteration,
            inner_iteration,
            residual,
            output,
        );
    }
}

/// Newmark kinematics from the acceleration-bearing correction. None of these
/// quantities subtracts two absolute positions or divides that difference by h².
fn structural_increment(
    method: &OperatorGeneralizedAlpha,
    v0: f64,
    a0: f64,
    correction: f64,
) -> (f64, f64, f64) {
    let a_new = correction / (method.beta * method.h * method.h);
    let dq = method.h.mul_add(
        v0,
        (method.h * method.h * (0.5 - method.beta)).mul_add(a0, correction),
    );
    let dv = method.h * (method.gamma.mul_add(a_new, (1.0 - method.gamma) * a0));
    (dq, dv, a_new)
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

/// Coefficients of the first-order step's effective Jacobian
/// `mass*M + tangent*Dr(t, u)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstOrderOperatorWeights {
    /// `alpha_m / (gamma*h)`.
    pub mass: f64,
    /// `alpha_f`.
    pub tangent: f64,
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

    /// Approximately solve the effective Jacobian system for `residual`.
    ///
    /// `t` and `u` are exactly the stage time and state used by `tangent_apply`;
    /// `weights` give the effective mass and tangent coefficients. The logical
    /// Newton/Krylov indices permit variable right preconditioning. Overwrite
    /// every output component with a finite value; invalid output refuses the
    /// solve before publication. The default is identity. Implementations must
    /// bound their own work and be deterministic for replay. The separately
    /// supplied adjoint transpose-solve preconditioner is unaffected.
    #[allow(clippy::too_many_arguments)] // stage, operator weights, and both logical iterations
    fn preconditioner_apply(
        &self,
        _t: f64,
        _u: &[f64],
        _weights: FirstOrderOperatorWeights,
        _outer_iteration: usize,
        _inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.copy_from_slice(residual);
    }
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
        self.step_controlled(state, problem, forcing, &mut || false)
    }

    /// Advance with cooperative cancellation before setup, before/after
    /// Newton residual, Jacobian and preconditioner callbacks, within inner
    /// Krylov and backtracking work, and before publication. Individual
    /// callbacks must bound their work. Cancellation leaves all state and
    /// history unchanged.
    pub fn step_controlled<P, Cancel>(
        &self,
        state: &mut FirstOrderState,
        problem: &P,
        forcing: &[f64],
        cancelled: &mut Cancel,
    ) -> Result<ImplicitStepTelemetry, TimeSolveError>
    where
        P: FirstOrderProblem + ?Sized,
        Cancel: FnMut() -> bool,
    {
        first_order_poll(cancelled)?;
        require_dimension("first-order problem", self.n, problem.dimension())?;
        require_dimension("first-order state", self.n, state.u.len())?;
        require_dimension("first-order rate", self.n, state.rate.len())?;
        require_dimension("forcing", self.n, forcing.len())?;
        let next_t = checked_time_advance(state.t, self.h)?;
        let t_eval = self.h.mul_add(self.alpha_f, state.t);
        if !(t_eval.is_finite() && t_eval > state.t) {
            return Err(TimeSolveError::InvalidStageTime {
                t_bits: state.t.to_bits(),
                h_bits: self.h.to_bits(),
                fraction_bits: self.alpha_f.to_bits(),
                stage_t_bits: t_eval.to_bits(),
            });
        }
        let next_steps = state
            .steps
            .checked_add(1)
            .ok_or(TimeSolveError::StepCounterOverflow { steps: state.steps })?;
        let residual = FirstOrderStepResidual {
            method: self,
            problem,
            t_eval,
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
        let setup = NewtonKrylovState::new(&residual, guess, self.solve.newton);
        first_order_poll(cancelled)?;
        let mut newton = setup.map_err(TimeSolveError::NewtonSetup)?;
        let mut report = newton.run(&residual, 0);
        for _ in 0..self.solve.max_newton_iterations {
            // Plateau is a report diagnosis, not a terminal Newton state;
            // the ordinary run continues after it, so this path must too.
            if report.converged
                || !matches!(
                    report.diagnosis,
                    Some(
                        fs_solver::NewtonStallDiagnosis::BudgetExhausted
                            | fs_solver::NewtonStallDiagnosis::Plateau
                    )
                )
            {
                break;
            }
            first_order_poll(cancelled)?;
            report = completed_newton(newton.run_controlled(&residual, 1, cancelled))?;
            first_order_poll(cancelled)?;
        }
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
        first_order_poll(cancelled)?;
        state.u = u_new;
        state.rate = rate_new;
        state.t = next_t;
        state.steps = next_steps;
        state.history.push(telemetry.clone());
        Ok(telemetry)
    }
}

fn first_order_poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), TimeSolveError> {
    if cancelled() {
        Err(TimeSolveError::Cancelled)
    } else {
        Ok(())
    }
}

struct FirstOrderStepResidual<'a, P: ?Sized> {
    method: &'a OperatorFirstOrderGeneralizedAlpha,
    problem: &'a P,
    t_eval: f64,
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
        let mut mass = vec![0.0; method.n];
        let mut internal = vec![0.0; method.n];
        self.problem.mass_apply(&rate_eval, &mut mass);
        self.problem
            .internal_force(self.t_eval, &u_eval, &mut internal);
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
        let mut tangent = vec![0.0; method.n];
        self.problem
            .tangent_apply(self.t_eval, &u_eval, &tangent_direction, &mut tangent);
        for i in 0..method.n {
            output[i] = mass_scale.mul_add(mass[i], tangent[i]);
        }
    }

    fn preconditioner_apply(
        &self,
        u_new: &[f64],
        outer_iteration: usize,
        inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        let method = self.method;
        let u_eval: Vec<f64> = u_new
            .iter()
            .zip(self.u0)
            .map(|(&u_new, &u0)| method.alpha_f.mul_add(u_new, (1.0 - method.alpha_f) * u0))
            .collect();
        let weights = FirstOrderOperatorWeights {
            mass: method.alpha_m / (method.gamma * method.h),
            tangent: method.alpha_f,
        };
        self.problem.preconditioner_apply(
            self.t_eval,
            &u_eval,
            weights,
            outer_iteration,
            inner_iteration,
            residual,
            output,
        );
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
