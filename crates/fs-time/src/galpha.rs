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
/// Fixed-configuration momentum jumps with consistent post-impulse acceleration.
#[path = "galpha/impulse.rs"]
pub mod impulse;
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
    /// A controlled linear action refused its shape or temporary allocation.
    SolverOperatorFailed(fs_solver::OperatorApplyError),
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
            Self::SolverOperatorFailed(error) => {
                write!(f, "generalized-alpha operator refused; step discarded: {error}")
            }
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
        SolverRunProgress::OperatorFailed(error) => {
            Err(TimeSolveError::SolverOperatorFailed(error))
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

    /// Overwrite with `R(q_base + fraction * dq, v) - forcing`.
    ///
    /// Split displacement coordinates retain a resolved stage increment even
    /// when adding it to the absolute displacement would round it away. The
    /// default reconstructs the same fused stage coordinate used by admission
    /// and derivative callbacks, then calls the comple…9058 tokens truncated…d, self.rate0[i], correction[i]);
            u_eval[i] = method.alpha_f.mul_add(increment, self.u0[i]);
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
        correction: &[f64],
        outer_iteration: usize,
        inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        let method = self.method;
        let u_eval: Vec<f64> = correction
            .iter()
            .enumerate()
            .map(|(i, &correction)| {
                let (increment, _) = first_order_increment(method, self.rate0[i], correction);
                method.alpha_f.mul_add(increment, self.u0[i])
            })
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

/// State increment and endpoint rate from the rate-bearing correction.
/// Recovering the rate never subtracts two rounded absolute states.
fn first_order_increment(
    method: &OperatorFirstOrderGeneralizedAlpha,
    rate0: f64,
    correction: f64,
) -> (f64, f64) {
    (
        (method.h * (1.0 - method.gamma)).mul_add(rate0, correction),
        correction / (method.gamma * method.h),
    )
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
