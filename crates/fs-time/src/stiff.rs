//! Stiff machinery: a second-order IMEX step (implicit stiff LINEAR
//! part via a prefactored LU, explicit nonlinearity) and an exponential
//! Euler for u′ = A·u + N(u) with SYMMETRIC A (eigenbasis via the landed
//! Jacobi; φ₁ evaluated as expm1(x)/x — the cancellation-free form
//! fs-math's expm1 exists for). Krylov φ-actions for large nonsymmetric
//! A are recorded follow-up (needs Arnoldi).

use fs_exec::reduce::Compensated;
use fs_ivl::expansion::fast_expansion_sum_zeroelim;
use fs_la::eigen::{MAX_EIGEN_WORK_ELEMENTS, admit_jacobi_eigh, jacobi_eigh};
use fs_la::factor::{Lu, lu};
use fs_math::det;
use fs_math::eft::two_sum;
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolveReport};
use std::fmt;

pub mod adjoint;

/// Prefactored operators for the IMEX-θ two-stage (ARS(2,2,2)-style)
/// scheme on u′ = L·u + N(u).
pub struct Imex2 {
    n: usize,
    h: f64,
    l_mat: Vec<f64>,
    solve_gamma: Lu,
    gamma: f64,
}

impl Imex2 {
    /// Build from the stiff linear operator (row-major n×n) and step h.
    ///
    /// # Panics
    /// If (I − γhL) is singular (h out of the scheme's range).
    #[must_use]
    pub fn new(l_mat: &[f64], n: usize, h: f64) -> Imex2 {
        let gamma = 1.0 - std::f64::consts::FRAC_1_SQRT_2; // ARS(2,2,2)
        let mut m = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let id = if i == j { 1.0 } else { 0.0 };
                m[i * n + j] = (-gamma * h).mul_add(l_mat[i * n + j], id);
            }
        }
        let solve_gamma = lu(&m, n).expect("(I - gamma h L) must be nonsingular");
        Imex2 {
            n,
            h,
            l_mat: l_mat.to_vec(),
            solve_gamma,
            gamma,
        }
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

/// One ARS(2,2,2) IMEX step: L treated implicitly (γ-diagonal, one
/// prefactored LU reused for both stages), N explicitly; second order
/// in BOTH parts and R(∞) = 0 on the stiff part. Updates `u` in place.
pub fn imex2_step<N: Fn(&[f64], &mut [f64])>(im: &Imex2, u: &mut [f64], nonlin: &N) {
    let (n, h, g) = (im.n, im.h, im.gamma);
    let mut nu = vec![0.0f64; n];
    nonlin(u, &mut nu);
    // Stage 1 (backward-Euler-γ on L, explicit N):
    // (I − γhL)·u₁ = u + γh·N(u).
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        rhs[i] = h.mul_add(g * nu[i], u[i]);
    }
    im.solve_gamma.solve(&mut rhs);
    let u1 = rhs;
    // Stage 2 (ARS(2,2,2), stiffly accurate — u⁺ IS the last stage):
    // (I − γhL)·u⁺ = u + h·[δ·N(u) + (1−δ)·N(u₁) + (1−γ)·L·u₁],
    // δ = 1 − 1/(2γ). The (δ, 1−δ) explicit weights — NOT trapezoidal
    // (½, ½), which drops the nonlinear part to first order — satisfy
    // (1−δ)γ = ½, the h²·N′N order condition.
    let delta = 1.0 - 1.0 / (2.0 * g);
    let mut nu1 = vec![0.0f64; n];
    nonlin(&u1, &mut nu1);
    let mut lu1 = vec![0.0f64; n];
    matvec(&im.l_mat, n, &u1, &mut lu1);
    let mut rhs2 = vec![0.0f64; n];
    for i in 0..n {
        let nbar = delta.mul_add(nu[i], (1.0 - delta) * nu1[i]);
        rhs2[i] = u[i] + h * ((1.0 - g) * lu1[i] + nbar);
    }
    im.solve_gamma.solve(&mut rhs2);
    u.copy_from_slice(&rhs2);
}

/// Krylov policy for each matrix-free IMEX stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImexSolveConfig {
    /// True relative-residual target.
    pub tolerance: f64,
    /// FGMRES restart length.
    pub restart: usize,
    /// Maximum restart cycles per stage.
    pub max_cycles: usize,
}

impl Default for ImexSolveConfig {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-11,
            restart: 24,
            max_cycles: 8,
        }
    }
}

/// Identity flexible preconditioner for small fixtures and unpreconditioned
/// operator probes. Field-scale callers should supply a problem-specific
/// `FlexiblePreconditioner` instead.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdentityPreconditioner;

impl FlexiblePreconditioner for IdentityPreconditioner {
    fn apply(&self, _logical_iteration: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

/// Complete iteration receipt for one accepted operator-backed IMEX step.
#[derive(Debug, Clone)]
pub struct ImexStepTelemetry {
    /// Zero-based accepted step index.
    pub step: usize,
    /// Time at the beginning of the step.
    pub t_start: f64,
    /// Step size.
    pub h: f64,
    /// First-stage FGMRES report.
    pub stage_one: SolveReport,
    /// Second-stage FGMRES report.
    pub stage_two: SolveReport,
}

impl PartialEq for ImexStepTelemetry {
    fn eq(&self, other: &Self) -> bool {
        self.step == other.step
            && self.t_start.to_bits() == other.t_start.to_bits()
            && self.h.to_bits() == other.h.to_bits()
            && solve_report_bits_equal(&self.stage_one, &other.stage_one)
            && solve_report_bits_equal(&self.stage_two, &other.stage_two)
    }
}

/// IMEX stage named by a refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImexStage {
    /// Solve of `(I - gamma h L) u1`.
    One,
    /// Solve of `(I - gamma h L) u_next`.
    Two,
}

/// Typed refusal from an operator-backed IMEX step.
#[derive(Debug, Clone)]
pub enum ImexSolveError {
    /// Cancellation was observed before publishing a complete step.
    Cancelled,
    /// State/operator dimensions disagree.
    Dimension {
        /// Required dimension.
        expected: usize,
        /// Supplied dimension.
        actual: usize,
    },
    /// The configured positive step cannot advance the supplied clock to a
    /// distinct finite time.
    InvalidTimeAdvance {
        /// Exact source-time bits.
        t_bits: u64,
        /// Exact configured-step bits.
        h_bits: u64,
        /// Exact rounded-sum bits.
        next_t_bits: u64,
    },
    /// The accepted-step counter cannot represent another commit.
    StepCounterOverflow {
        /// Exact current counter value.
        steps: usize,
    },
    /// The explicit nonlinearity produced NaN or infinity.
    NonFiniteNonlinearity {
        /// Stage whose explicit evaluation failed.
        stage: ImexStage,
        /// Entry index.
        index: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// One shifted linear solve exhausted its budget or broke down.
    NotConverged {
        /// Failed stage.
        stage: ImexStage,
        /// True-residual solve report.
        report: SolveReport,
    },
}

impl fmt::Display for ImexSolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cancelled => write!(f, "IMEX step cancelled before publication"),
            Self::Dimension { expected, actual } => write!(
                f,
                "IMEX operator/state dimension {actual} differs from method dimension {expected}"
            ),
            Self::InvalidTimeAdvance {
                t_bits,
                h_bits,
                next_t_bits,
            } => write!(
                f,
                "IMEX clock cannot advance finitely (t bits 0x{t_bits:016x}, h bits \
                 0x{h_bits:016x}, next bits 0x{next_t_bits:016x})"
            ),
            Self::StepCounterOverflow { steps } => {
                write!(f, "IMEX accepted-step counter cannot advance past {steps}")
            }
            Self::NonFiniteNonlinearity { stage, index, bits } => write!(
                f,
                "IMEX stage {stage:?} nonlinearity entry {index} is non-finite (bits 0x{bits:016x})"
            ),
            Self::NotConverged { stage, report } => write!(
                f,
                "IMEX stage {stage:?} did not converge after {} Krylov iterations: {:?}",
                report.iters, report.diagnosis
            ),
        }
    }
}

impl core::error::Error for ImexSolveError {}

/// Plain-data operator-backed IMEX trajectory checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct ImexState {
    /// Current time.
    pub t: f64,
    /// Current solution.
    pub u: Vec<f64>,
    /// Accepted steps.
    pub steps: usize,
    /// Complete accepted-step telemetry.
    pub history: Vec<ImexStepTelemetry>,
}

impl ImexState {
    /// Construct a checkpoint at one exact time.
    #[must_use]
    pub fn new(t: f64, u: &[f64]) -> Self {
        Self {
            t,
            u: u.to_vec(),
            steps: 0,
            history: Vec::new(),
        }
    }
}

/// ARS(2,2,2) over an arbitrary `LinearOp`, with injected flexible
/// preconditioning for both shifted stage systems.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OperatorImex2 {
    n: usize,
    h: f64,
    gamma: f64,
    solve: ImexSolveConfig,
}

struct ImexStages {
    one: Vec<f64>,
    next: Vec<f64>,
    report_one: SolveReport,
    report_two: SolveReport,
}

fn imex_poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), ImexSolveError> {
    if cancelled() {
        Err(ImexSolveError::Cancelled)
    } else {
        Ok(())
    }
}

impl OperatorImex2 {
    /// Configure an operator-backed IMEX method.
    #[must_use]
    pub fn new(n: usize, h: f64, solve: ImexSolveConfig) -> Self {
        assert!(n > 0, "IMEX dimension must be positive");
        assert!(
            h.is_finite() && h > 0.0,
            "IMEX h must be positive and finite"
        );
        assert!(
            solve.tolerance.is_finite() && solve.tolerance > 0.0,
            "IMEX tolerance must be positive and finite"
        );
        assert!(solve.restart > 0, "IMEX restart must be positive");
        assert!(solve.max_cycles > 0, "IMEX cycle budget must be positive");
        Self {
            n,
            h,
            gamma: 1.0 - std::f64::consts::FRAC_1_SQRT_2,
            solve,
        }
    }

    /// Advance one step. State changes only after both true-residual solves
    /// converge, so a failed attempt is transaction-like.
    pub fn step<L, N, P>(
        &self,
        state: &mut ImexState,
        linear: &L,
        preconditioner: &P,
        nonlin: &N,
    ) -> Result<ImexStepTelemetry, ImexSolveError>
    where
        L: LinearOp + ?Sized,
        N: Fn(&[f64], &mut [f64]),
        P: FlexiblePreconditioner,
    {
        self.step_controlled(state, linear, preconditioner, nonlin, &mut || false)
    }

    /// Cancellable ARS(2,2,2) step. Polls around model evaluations, between
    /// FGMRES restart cycles, and immediately before publication. One cycle
    /// contains at most the configured `restart` Arnoldi iterations; callbacks
    /// must themselves bound their work. Cancellation leaves the entire state,
    /// including its clock and history, unchanged, so the attempt is retryable.
    pub fn step_controlled<L, N, P, Cancel>(
        &self,
        state: &mut ImexState,
        linear: &L,
        preconditioner: &P,
        nonlin: &N,
        cancelled: &mut Cancel,
    ) -> Result<ImexStepTelemetry, ImexSolveError>
    where
        L: LinearOp + ?Sized,
        N: Fn(&[f64], &mut [f64]),
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        imex_poll(cancelled)?;
        if linear.n() != self.n {
            return Err(ImexSolveError::Dimension {
                expected: self.n,
                actual: linear.n(),
            });
        }
        if state.u.len() != self.n {
            return Err(ImexSolveError::Dimension {
                expected: self.n,
                actual: state.u.len(),
            });
        }
        let next_t = state.t + self.h;
        if !(state.t.is_finite() && next_t.is_finite() && next_t > state.t) {
            return Err(ImexSolveError::InvalidTimeAdvance {
                t_bits: state.t.to_bits(),
                h_bits: self.h.to_bits(),
                next_t_bits: next_t.to_bits(),
            });
        }
        let next_steps = state
            .steps
            .checked_add(1)
            .ok_or(ImexSolveError::StepCounterOverflow { steps: state.steps })?;
        let stages = self.stages(&state.u, linear, preconditioner, nonlin, cancelled)?;
        let telemetry = ImexStepTelemetry {
            step: state.steps,
            t_start: state.t,
            h: self.h,
            stage_one: stages.report_one,
            stage_two: stages.report_two,
        };
        imex_poll(cancelled)?;
        state.u = stages.next;
        state.t = next_t;
        state.steps = next_steps;
        state.history.push(telemetry.clone());
        Ok(telemetry)
    }

    fn solve_stage<L, P, Cancel>(
        &self,
        operator: &L,
        preconditioner: &P,
        rhs: &[f64],
        stage: ImexStage,
        cancelled: &mut Cancel,
    ) -> Result<(Vec<f64>, SolveReport), ImexSolveError>
    where
        L: LinearOp,
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        imex_poll(cancelled)?;
        let mut solve = FgmresState::new(rhs, self.solve.restart);
        for cycle in 0..self.solve.max_cycles {
            imex_poll(cancelled)?;
            let report = solve.run(operator, preconditioner, rhs, self.solve.tolerance, 1);
            imex_poll(cancelled)?;
            if report.converged {
                return Ok((solve.x, report));
            }
            if report.diagnosis == Some(fs_solver::StallDiagnosis::Breakdown)
                || cycle + 1 == self.solve.max_cycles
            {
                return Err(ImexSolveError::NotConverged { stage, report });
            }
        }
        unreachable!("constructor requires a positive cycle budget")
    }

    fn stages<L, N, P, Cancel>(
        &self,
        initial: &[f64],
        linear: &L,
        preconditioner: &P,
        nonlin: &N,
        cancelled: &mut Cancel,
    ) -> Result<ImexStages, ImexSolveError>
    where
        L: LinearOp + ?Sized,
        N: Fn(&[f64], &mut [f64]),
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        imex_poll(cancelled)?;
        let shifted = ShiftedLinearOp {
            linear,
            shift: -self.gamma * self.h,
        };
        let nu = evaluate_nonlinearity(nonlin, initial, ImexStage::One)?;
        imex_poll(cancelled)?;
        let mut rhs_one = vec![0.0; self.n];
        for i in 0..self.n {
            rhs_one[i] = self.h.mul_add(self.gamma * nu[i], initial[i]);
        }
        let (u_one, report_one) = self.solve_stage(
            &shifted,
            preconditioner,
            &rhs_one,
            ImexStage::One,
            cancelled,
        )?;
        let nu_one = evaluate_nonlinearity(nonlin, &u_one, ImexStage::Two)?;
        imex_poll(cancelled)?;
        let mut linear_u_one = vec![0.0; self.n];
        linear.apply(&u_one, &mut linear_u_one);
        imex_poll(cancelled)?;
        let delta = 1.0 - 1.0 / (2.0 * self.gamma);
        let mut rhs_two = vec![0.0; self.n];
        for i in 0..self.n {
            let explicit = delta.mul_add(nu[i], (1.0 - delta) * nu_one[i]);
            rhs_two[i] = initial[i] + self.h * ((1.0 - self.gamma) * linear_u_one[i] + explicit);
        }
        let (next, report_two) = self.solve_stage(
            &shifted,
            preconditioner,
            &rhs_two,
            ImexStage::Two,
            cancelled,
        )?;
        Ok(ImexStages {
            one: u_one,
            next,
            report_one,
            report_two,
        })
    }
}

struct ShiftedLinearOp<'a, L: ?Sized> {
    linear: &'a L,
    shift: f64,
}

impl<L: LinearOp + ?Sized> LinearOp for ShiftedLinearOp<'_, L> {
    fn n(&self) -> usize {
        self.linear.n()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) {
        self.linear.apply(input, output);
        for (value, input_value) in output.iter_mut().zip(input) {
            *value = self.shift.mul_add(*value, *input_value);
        }
    }

    fn apply_transpose(&self, input: &[f64], output: &mut [f64]) {
        self.linear.apply_transpose(input, output);
        for (value, input_value) in output.iter_mut().zip(input) {
            *value = self.shift.mul_add(*value, *input_value);
        }
    }
}

fn evaluate_nonlinearity<N: Fn(&[f64], &mut [f64])>(
    nonlin: &N,
    state: &[f64],
    stage: ImexStage,
) -> Result<Vec<f64>, ImexSolveError> {
    let mut output = vec![0.0; state.len()];
    nonlin(state, &mut output);
    for (index, value) in output.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(ImexSolveError::NonFiniteNonlinearity {
                stage,
                index,
                bits: value.to_bits(),
            });
        }
    }
    Ok(output)
}

fn solve_report_bits_equal(left: &SolveReport, right: &SolveReport) -> bool {
    left.iters == right.iters
        && left.rel_residual.to_bits() == right.rel_residual.to_bits()
        && left.converged == right.converged
        && left.diagnosis == right.diagnosis
        && left.history.len() == right.history.len()
        && left
            .history
            .iter()
            .zip(&right.history)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

/// Exponential Euler for u′ = A·u + N(u), SYMMETRIC A:
/// u⁺ = e^{hA}·u + h·φ₁(hA)·N(u), computed in A's eigenbasis with
/// φ₁(x) = expm1(x)/x (exact limit 1 at x = 0). Exact for constant N,
/// up to floating-point and eigensolve error; first order for nonlinear N.
/// Preparation scales the matrix and checks the resulting eigenpairs. It is
/// a bounded dense kernel, with no internal cancellation point or certificate
/// of the physical model. Use [`Self::workspace`] and [`Self::step_checked`]
/// to reuse storage and refuse invalid steps without changing the state.
pub struct ExpEuler {
    n: usize,
    h: f64,
    /// Eigenvectors (columns) of A.
    vecs: Vec<f64>,
    /// e^{hλ}, or e^{hλ/4} when direct retention is subnormal/zero.
    exp_h: Vec<f64>,
    /// Four successive factors retain finite tails of large initial states.
    split_exp: Vec<bool>,
    /// h·φ₁(hλ) per eigenvalue.
    hphi1: Vec<f64>,
    /// Exact constant modes of explicitly reciprocal diffusion components.
    constant_modes: Vec<ConstantDiffusionMode>,
    /// Component nodes in the same order as `constant_modes`.
    component_nodes: Vec<usize>,
}

struct ConstantDiffusionMode {
    column: usize,
    start: usize,
    end: usize,
    normalization: f64,
}

/// A refusal from spectral preparation or exponential stepping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExpEulerError {
    /// Invalid step, matrix, or dense workspace admission.
    InvalidInput(&'static str),
    /// A matrix, state, or reusable workspace has the wrong length.
    Dimension {
        /// Required length.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// Opposite matrix entries differ by more than roundoff.
    NonSymmetric {
        /// Row of the upper-triangle entry.
        row: usize,
        /// Column of the upper-triangle entry.
        column: usize,
    },
    /// The eigenpairs or scalar response could not be resolved in binary64.
    Spectrum(&'static str),
    /// An input state entry is nonfinite.
    NonFiniteState {
        /// Entry index.
        index: usize,
    },
    /// The callback left an entry unwritten or supplied a nonfinite source.
    NonFiniteSource {
        /// Entry index.
        index: usize,
    },
    /// A modal projection or proposed endpoint is unrepresentable.
    Arithmetic,
}

impl fmt::Display for ExpEulerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-TIME-EXP-EULER: {self:?}")
    }
}

impl std::error::Error for ExpEulerError {}

/// Reusable scratch storage for one exponential step. It owns no accepted
/// state and can be reused with any [`ExpEuler`] of the same dimension.
pub struct ExpEulerWorkspace {
    source: Vec<f64>,
    modal: Vec<f64>,
    next: Vec<f64>,
}

impl ExpEuler {
    /// Build from symmetric A (row-major n×n) and step h.
    ///
    /// # Panics
    /// Panics on a refusal from [`Self::try_new`].
    #[must_use]
    pub fn new(a: &[f64], n: usize, h: f64) -> ExpEuler {
        Self::try_new(a, n, h).expect("valid symmetric exponential-Euler problem")
    }

    /// Prepare a finite symmetric operator and positive finite duration.
    /// Opposite entries may differ by at most 64 eps relative to that pair;
    /// their mean then defines the symmetric operator. The scaled Jacobi
    /// eigenpairs must pass residual and orthogonality checks at roundoff
    /// scale. Unresolved dynamic range and unrepresentable growth refuse.
    /// Exactly reciprocal diffusion matrices with nonnegative off-diagonal
    /// entries and exactly zero row sums use [`Self::try_new_diffusion`]'s
    /// conserved modes. Exact expansion sums distinguish a true nullspace
    /// from a small physical growth or decay rate; no eigenvalue is clipped.
    ///
    /// # Errors
    /// Returns [`ExpEulerError`] before creating a usable stepper when the
    /// inputs, spectrum, or scalar propagator cannot be admitted.
    pub fn try_new(a: &[f64], n: usize, h: f64) -> Result<Self, ExpEulerError> {
        Self::try_new_symmetric(a, n, h, false)
    }

    fn try_new_symmetric(
        a: &[f64],
        n: usize,
        h: f64,
        require_decay: bool,
    ) -> Result<Self, ExpEulerError> {
        if !h.is_finite() || h <= 0.0 {
            return Err(ExpEulerError::InvalidInput(
                "require a positive finite step",
            ));
        }
        let admission = admit_jacobi_eigh(n)
            .map_err(|_| ExpEulerError::InvalidInput("dense eigensolve exceeds its work cap"))?;
        let square = admission.matrix_entries();
        if a.len() != square {
            return Err(ExpEulerError::Dimension {
                expected: square,
                actual: a.len(),
            });
        }
        if a.iter().any(|value| !value.is_finite()) {
            return Err(ExpEulerError::InvalidInput("require a finite matrix"));
        }
        if !require_decay && n > 0 && exact_diffusion_generator(a, n) {
            return Self::prepare_diffusion(a, n, h);
        }
        let maximum = a
            .iter()
            .fold(0.0_f64, |scale, value| scale.max(value.abs()));
        // Binary scaling changes no ordinary-range significand bits. Besides
        // avoiding diagonal-sum overflow inside Jacobi, it preserves the
        // existing operation tree for already well-scaled operators.
        let exponent_bits = maximum.to_bits() & 0x7ff0_0000_0000_0000;
        let scale = if exponent_bits == 0 {
            f64::MIN_POSITIVE
        } else {
            f64::from_bits(exponent_bits)
        };
        let mut scaled = vec![0.0; square];
        for i in 0..n {
            for j in i..n {
                let upper = a[i * n + j];
                let lower = a[j * n + i];
                let pair_scale = upper.abs().max(lower.abs());
                if pair_scale > 0.0
                    && (upper / pair_scale - lower / pair_scale).abs() > 64.0 * f64::EPSILON
                {
                    return Err(ExpEulerError::NonSymmetric { row: i, column: j });
                }
                let value = if scale == 0.0 {
                    0.0
                } else {
                    0.5 * (upper / scale + lower / scale)
                };
                if pair_scale > 0.0 && value == 0.0 {
                    return Err(ExpEulerError::Spectrum("matrix scaling underflows"));
                }
                scaled[i * n + j] = value;
                scaled[j * n + i] = value;
            }
        }
        let (vals, vecs) = jacobi_eigh(&scaled, n);
        let tolerance = 128.0 * f64::EPSILON * n.max(1) as f64;
        let norm = scaled.chunks(n.max(1)).fold(0.0_f64, |norm, row| {
            norm.max(row.iter().map(|value| value.abs()).sum())
        });
        for i in 0..n {
            if !vals[i].is_finite() {
                return Err(ExpEulerError::Spectrum("nonfinite eigenvalue"));
            }
            for j in 0..n {
                let mut av = 0.0;
                let mut gram = 0.0;
                for k in 0..n {
                    av = scaled[j * n + k].mul_add(vecs[k * n + i], av);
                    gram = vecs[k * n + j].mul_add(vecs[k * n + i], gram);
                }
                let residual = av - vals[i] * vecs[j * n + i];
                let identity = if i == j { 1.0 } else { 0.0 };
                if !residual.is_finite()
                    || residual.abs() > tolerance * norm
                    || !gram.is_finite()
                    || (gram - identity).abs() > tolerance
                {
                    return Err(ExpEulerError::Spectrum("unresolved eigenbasis"));
                }
            }
        }
        let mut exp_h = Vec::with_capacity(n);
        let mut split_exp = Vec::with_capacity(n);
        let mut hphi1 = Vec::with_capacity(n);
        for value in vals {
            let lambda = value * scale;
            if !lambda.is_finite() || (value != 0.0 && scale > 0.0 && lambda == 0.0) {
                return Err(ExpEulerError::Spectrum(
                    "eigenvalue scaling is unrepresentable",
                ));
            }
            if require_decay && lambda >= 0.0 {
                return Err(ExpEulerError::Spectrum(
                    "unresolved diffusion contrast mode",
                ));
            }
            let x = h * lambda;
            let retention = det::exp(x);
            // For stiff decay, dividing by lambda avoids h*(1/|h*lambda|):
            // the latter erases a finite forced equilibrium when h*lambda
            // overflows. Near zero use phi_1 to retain its exact limit.
            let source = if x.abs() < f64::EPSILON {
                h
            } else if x.abs() <= 1.0 {
                h * (det::expm1(x) / x)
            } else {
                det::expm1(x) / lambda
            };
            if !retention.is_finite() || !source.is_finite() || source <= 0.0 {
                return Err(ExpEulerError::Spectrum(
                    "scalar response is unrepresentable",
                ));
            }
            let split = retention < f64::MIN_POSITIVE;
            exp_h.push(if split { det::exp(0.25 * x) } else { retention });
            split_exp.push(split);
            hphi1.push(source);
        }
        Ok(Self {
            n,
            h,
            vecs,
            exp_h,
            split_exp,
            hphi1,
            constant_modes: Vec::new(),
            component_nodes: Vec::new(),
        })
    }

    /// Prepare reciprocal diffusion from a row-major symmetric rate table.
    ///
    /// The diagonal must be zero and off-diagonal rates must be finite,
    /// nonnegative and exactly reciprocal. They define
    /// `u_i' = sum_j rate_ij * (u_j - u_i) + N_i(u)` for equal-capacity
    /// states. Rates have units of inverse time; geometry, diffusivity and
    /// boundary conductance remain the caller's physical model. No transport
    /// across absent edges or unequal-capacity interpretation is inferred.
    ///
    /// Each connected component retains its constant mode analytically. Only
    /// its orthogonal contrasts enter the eigensolver, so eigenvalue roundoff
    /// cannot turn conserved material into growth or decay at long times.
    /// Compensated component sums also retain balanced sources and sinks.
    /// The existing workspace and allocation-free step path are reused.
    /// Constant additive forcing is exact up to eigensolve and roundoff error;
    /// nonlinear forcing remains first order. This bounded dense owner adds
    /// no positivity, sparse-field scaling or internal cancellation guarantee.
    ///
    /// # Errors
    /// Refuses invalid rates, durations or shapes, the aggregate dense work
    /// cap, and unresolved/nonrepresentable decaying contrast spectra.
    pub fn try_new_diffusion(rates: &[f64], n: usize, h: f64) -> Result<Self, ExpEulerError> {
        let square = admit_jacobi_eigh(n)
            .map_err(|_| ExpEulerError::InvalidInput("dense eigensolve exceeds its work cap"))?
            .matrix_entries();
        if rates.len() != square {
            return Err(ExpEulerError::Dimension {
                expected: square,
                actual: rates.len(),
            });
        }
        if (0..n).any(|i| rates[i * n + i] != 0.0) {
            return Err(ExpEulerError::InvalidInput(
                "exchange-rate diagonal must be zero",
            ));
        }
        Self::prepare_diffusion(rates, n, h)
    }

    // The diagonal is ignored here: the public rate-table constructor checks
    // zero entries, while the matrix constructor proves exact row closure.
    fn prepare_diffusion(rates: &[f64], n: usize, h: f64) -> Result<Self, ExpEulerError> {
        if n == 0 || !h.is_finite() || h <= 0.0 {
            return Err(ExpEulerError::InvalidInput(
                "require positive dimension and finite step",
            ));
        }
        let square = admit_jacobi_eigh(n)
            .map_err(|_| ExpEulerError::InvalidInput("dense eigensolve exceeds its work cap"))?
            .matrix_entries();
        // The complete output basis, contrast matrix and original eigensolver
        // coexist. Count their aggregate payload before constructing any graph.
        if square
            .checked_mul(6)
            .and_then(|value| {
                n.checked_mul(12)
                    .and_then(|linear| value.checked_add(linear))
            })
            .is_none_or(|value| value > MAX_EIGEN_WORK_ELEMENTS)
        {
            return Err(ExpEulerError::InvalidInput(
                "diffusion eigensolve exceeds its work cap",
            ));
        }
        for i in 0..n {
            for j in i + 1..n {
                let rate = rates[i * n + j];
                if !rate.is_finite() || rate < 0.0 {
                    return Err(ExpEulerError::InvalidInput(
                        "require finite nonnegative exchange rates",
                    ));
                }
                if rate != rates[j * n + i] {
                    return Err(ExpEulerError::NonSymmetric { row: i, column: j });
                }
            }
        }
        let mut result = Self {
            n,
            h,
            vecs: vec![0.0; square],
            exp_h: vec![1.0; n],
            split_exp: vec![false; n],
            hphi1: vec![h; n],
            constant_modes: Vec::new(),
            component_nodes: Vec::with_capacity(n),
        };
        let mut visited = vec![false; n];
        let mut column = 0;
        for root in 0..n {
            if visited[root] {
                continue;
            }
            let start = result.component_nodes.len();
            result.component_nodes.push(root);
            visited[root] = true;
            let mut cursor = start;
            while cursor < result.component_nodes.len() {
                let node = result.component_nodes[cursor];
                for next in 0..n {
                    if !visited[next] && rates[node * n + next] > 0.0 {
                        visited[next] = true;
                        result.component_nodes.push(next);
                    }
                }
                cursor += 1;
            }
            let end = result.component_nodes.len();
            result.component_nodes[start..end].sort_unstable();
            let nodes = &result.component_nodes[start..end];
            let count = nodes.len();
            let root_count = det::sqrt(count as f64);
            let normalization = 1.0 / root_count;
            result.constant_modes.push(ConstantDiffusionMode {
                column,
                start,
                end,
                normalization,
            });
            for &node in nodes {
                result.vecs[node * n + column] = normalization;
            }
            column += 1;
            if count == 1 {
                continue;
            }
            let mut diagonal = Vec::with_capacity(count);
            for &node in nodes {
                let sum = nodes
                    .iter()
                    .filter(|&&other| other != node)
                    .fold(Compensated::zero(), |sum, &other| {
                        sum.accumulate(rates[node * n + other])
                    });
                diagonal.push(-sum.value());
            }
            // H maps the first coordinate to minus the normalized constant.
            // With A*1 = 0 known from pair exchange, the contrast block of
            // H A H is A_ij - b(A_i0 + A_0j) + b² A_00, b=1/(sqrt(m)+1).
            // The known zero row/column never enters a numerical eigensolve.
            let b = 1.0 / (root_count + 1.0);
            let dimension = count - 1;
            let mut contrast = vec![0.0; dimension * dimension];
            for i in 1..count {
                for j in i..count {
                    let entry = if i == j {
                        diagonal[i]
                    } else {
                        rates[nodes[i] * n + nodes[j]]
                    };
                    let root_edges =
                        rates[nodes[i] * n + nodes[0]] + rates[nodes[0] * n + nodes[j]];
                    let value = (b * b).mul_add(diagonal[0], (-b).mul_add(root_edges, entry));
                    contrast[(i - 1) * dimension + j - 1] = value;
                    contrast[(j - 1) * dimension + i - 1] = value;
                }
            }
            let reduced = Self::try_new_symmetric(&contrast, dimension, h, true)?;
            for mode in 0..dimension {
                let sum = (0..dimension).fold(Compensated::zero(), |sum, row| {
                    sum.accumulate(reduced.vecs[row * dimension + mode])
                });
                let first = -normalization * sum.value();
                result.vecs[nodes[0] * n + column] = first;
                for row in 1..count {
                    result.vecs[nodes[row] * n + column] =
                        b.mul_add(first, reduced.vecs[(row - 1) * dimension + mode]);
                }
                result.exp_h[column] = reduced.exp_h[mode];
                result.split_exp[column] = reduced.split_exp[mode];
                result.hphi1[column] = reduced.hphi1[mode];
                column += 1;
            }
        }
        Ok(result)
    }

    /// The step size.
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// Allocate scratch once for repeated allocation-free checked steps.
    #[must_use]
    pub fn workspace(&self) -> ExpEulerWorkspace {
        ExpEulerWorkspace {
            source: vec![0.0; self.n],
            modal: vec![0.0; self.n],
            next: vec![0.0; self.n],
        }
    }

    /// One exponential-Euler step, allocating temporary workspace. The callback
    /// must write every source entry, including explicit zeros.
    ///
    /// # Panics
    /// Panics on a refusal from [`Self::step_checked`], before publishing state.
    pub fn step<N: Fn(&[f64], &mut [f64])>(&self, u: &mut [f64], nonlin: &N) {
        self.step_checked(u, nonlin, &mut self.workspace())
            .expect("valid exponential-Euler state and source");
    }

    /// Advance with reusable storage. The callback reads the original state
    /// and must write every source entry. No allocation occurs in this method;
    /// callback side effects and allocation remain the caller's responsibility.
    ///
    /// # Errors
    /// Invalid shapes, nonfinite/unwritten values and arithmetic overflow
    /// return an error without changing any entry of `u`. Workspace contents
    /// are scratch and may change on a refusal.
    pub fn step_checked<N: Fn(&[f64], &mut [f64])>(
        &self,
        u: &mut [f64],
        nonlin: &N,
        workspace: &mut ExpEulerWorkspace,
    ) -> Result<(), ExpEulerError> {
        let n = self.n;
        for actual in [u.len(), workspace.source.len()] {
            if actual != n {
                return Err(ExpEulerError::Dimension {
                    expected: n,
                    actual,
                });
            }
        }
        for (index, value) in u.iter().enumerate() {
            if !value.is_finite() {
                return Err(ExpEulerError::NonFiniteState { index });
            }
        }
        workspace.source.fill(f64::NAN);
        nonlin(u, &mut workspace.source);
        for (index, value) in workspace.source.iter().enumerate() {
            if !value.is_finite() {
                return Err(ExpEulerError::NonFiniteSource { index });
            }
        }
        // Transform to eigenbasis: û = Vᵀu, n̂ = VᵀN.
        let mut constants = self.constant_modes.iter().peekable();
        for i in 0..n {
            if constants.peek().is_some_and(|mode| mode.column == i) {
                let mode = constants.next().expect("constant mode was matched");
                let nodes = &self.component_nodes[mode.start..mode.end];
                // Exact zero-eliminating sums preserve cancellation across
                // arbitrarily mixed finite magnitudes. A compensated scalar
                // sum can still invent a small source which a long h amplifies.
                // Endpoint scratch is unused until all modal values exist.
                let state = component_sum(u, nodes, &mut workspace.next)?;
                let source = component_sum(&workspace.source, nodes, &mut workspace.next)?;
                workspace.modal[i] = mode.normalization * self.h.mul_add(source, state);
                if !workspace.modal[i].is_finite() {
                    return Err(ExpEulerError::Arithmetic);
                }
                continue;
            }
            let (mut au, mut an) = (0.0f64, 0.0f64);
            for j in 0..n {
                au = self.vecs[j * n + i].mul_add(u[j], au);
                an = self.vecs[j * n + i].mul_add(workspace.source[j], an);
            }
            if !au.is_finite() || !an.is_finite() {
                return Err(ExpEulerError::Arithmetic);
            }
            workspace.modal[i] = if self.split_exp[i] {
                // A vanished e^(h lambda) can still leave a finite tail when
                // u is large. Four factors stay normal throughout that range.
                let factor = self.exp_h[i];
                let tail = ((au * factor) * factor) * factor;
                factor.mul_add(tail, self.hphi1[i] * an)
            } else {
                self.exp_h[i].mul_add(au, self.hphi1[i] * an)
            };
            if !workspace.modal[i].is_finite() {
                return Err(ExpEulerError::Arithmetic);
            }
        }
        // Transform into scratch; publish only a complete finite endpoint.
        for j in 0..n {
            let mut acc = 0.0f64;
            for (i, &uhi) in workspace.modal.iter().enumerate() {
                acc = self.vecs[j * n + i].mul_add(uhi, acc);
            }
            if !acc.is_finite() {
                return Err(ExpEulerError::Arithmetic);
            }
            workspace.next[j] = acc;
        }
        u.copy_from_slice(&workspace.next);
        Ok(())
    }
}

// Grow a nonoverlapping expansion in existing endpoint scratch, with the
// shared error-free TwoSum primitive. Each input adds at most one component;
// writes never overtake unread components, so nodes.len() slots suffice.
// Exact balance becomes the empty expansion before any binary64 estimate.
// Refuse overflow of an intermediate expansion instead of publishing a false
// conserved amount. The caller owns finite inputs and scratch of length >= n.
fn component_sum(
    values: &[f64],
    nodes: &[usize],
    scratch: &mut [f64],
) -> Result<f64, ExpEulerError> {
    let mut length = 0;
    for &node in nodes {
        let mut carry = values[node];
        let mut output = 0;
        for input in 0..length {
            let (sum, error) = two_sum(carry, scratch[input]);
            if !sum.is_finite() || !error.is_finite() {
                return Err(ExpEulerError::Arithmetic);
            }
            if error != 0.0 {
                scratch[output] = error;
                output += 1;
            }
            carry = sum;
        }
        if carry != 0.0 {
            scratch[output] = carry;
            output += 1;
        }
        length = output;
    }
    let sum = fs_ivl::expansion::estimate(&scratch[..length]);
    if !sum.is_finite() {
        return Err(ExpEulerError::Arithmetic);
    }
    Ok(sum)
}

// Recognize only an exact reciprocal diffusion generator, not a near-zero
// eigenvalue. Expansion arithmetic preserves the distinction between a closed
// row and a physical loss smaller than a compensated reduction can resolve.
fn exact_diffusion_generator(a: &[f64], n: usize) -> bool {
    for i in 0..n {
        if a[i * n + i] > 0.0 {
            return false;
        }
        for j in i + 1..n {
            if a[i * n + j] < 0.0 || a[i * n + j] != a[j * n + i] {
                return false;
            }
        }
    }
    for i in 0..n {
        let diagonal = a[i * n + i];
        let mut sum = Compensated::zero().accumulate(diagonal);
        for j in 0..n {
            if i != j {
                sum = sum.accumulate(a[i * n + j]);
            }
        }
        if sum.value() != 0.0 {
            return false;
        }
        // Negative diagonal first, followed only by positive entries. The
        // finite near-zero fast sum above excludes overflow of this monotone
        // exact accumulation; no unrepresentable expansion is attempted.
        let mut exact = if diagonal == 0.0 {
            Vec::new()
        } else {
            vec![diagonal]
        };
        for j in 0..n {
            if i != j && a[i * n + j] != 0.0 {
                exact = fast_expansion_sum_zeroelim(&exact, &[a[i * n + j]]);
            }
        }
        if exact.iter().any(|&value| value != 0.0) {
            return false;
        }
    }
    true
}
