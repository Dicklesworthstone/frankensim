//! Stiff machinery: a second-order IMEX step (implicit stiff LINEAR
//! part via a prefactored LU, explicit nonlinearity) and an exponential
//! Euler for u′ = A·u + N(u) with SYMMETRIC A (eigenbasis via the landed
//! Jacobi; φ₁ evaluated as expm1(x)/x — the cancellation-free form
//! fs-math's expm1 exists for). Krylov φ-actions for large nonsymmetric
//! A are recorded follow-up (needs Arnoldi).

use fs_la::eigen::jacobi_eigh;
use fs_la::factor::{Lu, lu};
use fs_math::det;
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolveReport};
use std::fmt;

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
    /// State/operator dimensions disagree.
    Dimension {
        /// Required dimension.
        expected: usize,
        /// Supplied dimension.
        actual: usize,
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
            Self::Dimension { expected, actual } => write!(
                f,
                "IMEX operator/state dimension {actual} differs from method dimension {expected}"
            ),
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
    #[allow(clippy::too_many_lines)] // Both tableau stages form one atomic transaction.
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
        let shifted = ShiftedLinearOp {
            linear,
            shift: -self.gamma * self.h,
        };
        let nu = evaluate_nonlinearity(nonlin, &state.u, ImexStage::One)?;
        let mut rhs_one = vec![0.0; self.n];
        for i in 0..self.n {
            rhs_one[i] = self.h.mul_add(self.gamma * nu[i], state.u[i]);
        }
        let mut stage_one = FgmresState::new(&rhs_one, self.solve.restart);
        let report_one = stage_one.run(
            &shifted,
            preconditioner,
            &rhs_one,
            self.solve.tolerance,
            self.solve.max_cycles,
        );
        if !report_one.converged {
            return Err(ImexSolveError::NotConverged {
                stage: ImexStage::One,
                report: report_one,
            });
        }
        let u_one = stage_one.x;
        let nu_one = evaluate_nonlinearity(nonlin, &u_one, ImexStage::Two)?;
        let mut linear_u_one = vec![0.0; self.n];
        linear.apply(&u_one, &mut linear_u_one);
        let delta = 1.0 - 1.0 / (2.0 * self.gamma);
        let mut rhs_two = vec![0.0; self.n];
        for i in 0..self.n {
            let explicit = delta.mul_add(nu[i], (1.0 - delta) * nu_one[i]);
            rhs_two[i] = state.u[i] + self.h * ((1.0 - self.gamma) * linear_u_one[i] + explicit);
        }
        let mut stage_two = FgmresState::new(&rhs_two, self.solve.restart);
        let report_two = stage_two.run(
            &shifted,
            preconditioner,
            &rhs_two,
            self.solve.tolerance,
            self.solve.max_cycles,
        );
        if !report_two.converged {
            return Err(ImexSolveError::NotConverged {
                stage: ImexStage::Two,
                report: report_two,
            });
        }
        let telemetry = ImexStepTelemetry {
            step: state.steps,
            t_start: state.t,
            h: self.h,
            stage_one: report_one,
            stage_two: report_two,
        };
        state.u = stage_two.x;
        state.t += self.h;
        state.steps += 1;
        state.history.push(telemetry.clone());
        Ok(telemetry)
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
/// φ₁(x) = expm1(x)/x (exact limit 1 at x = 0). EXACT for N ≡ 0.
pub struct ExpEuler {
    n: usize,
    h: f64,
    /// Eigenvectors (columns) of A.
    vecs: Vec<f64>,
    /// e^{hλ} per eigenvalue.
    exp_h: Vec<f64>,
    /// h·φ₁(hλ) per eigenvalue.
    hphi1: Vec<f64>,
}

impl ExpEuler {
    /// Build from symmetric A (row-major n×n) and step h.
    #[must_use]
    pub fn new(a: &[f64], n: usize, h: f64) -> ExpEuler {
        let (vals, vecs) = jacobi_eigh(a, n);
        let exp_h: Vec<f64> = vals.iter().map(|&l| det::exp(h * l)).collect();
        let hphi1: Vec<f64> = vals
            .iter()
            .map(|&l| {
                let x = h * l;
                if x.abs() < 1e-300 {
                    h
                } else {
                    h * (det::expm1(x) / x)
                }
            })
            .collect();
        ExpEuler {
            n,
            h,
            vecs,
            exp_h,
            hphi1,
        }
    }

    /// The step size.
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// One exponential-Euler step (u updated in place).
    pub fn step<N: Fn(&[f64], &mut [f64])>(&self, u: &mut [f64], nonlin: &N) {
        let n = self.n;
        let mut nu = vec![0.0f64; n];
        nonlin(u, &mut nu);
        // Transform to eigenbasis: û = Vᵀu, n̂ = VᵀN.
        let mut uh = vec![0.0f64; n];
        let mut nh = vec![0.0f64; n];
        for i in 0..n {
            let (mut au, mut an) = (0.0f64, 0.0f64);
            for j in 0..n {
                au = self.vecs[j * n + i].mul_add(u[j], au);
                an = self.vecs[j * n + i].mul_add(nu[j], an);
            }
            uh[i] = au;
            nh[i] = an;
        }
        // Apply the scalar filters and transform back.
        for (i, uhi) in uh.iter_mut().enumerate() {
            *uhi = self.exp_h[i].mul_add(*uhi, self.hphi1[i] * nh[i]);
        }
        for (j, uj) in u.iter_mut().enumerate() {
            let mut acc = 0.0f64;
            for (i, &uhi) in uh.iter().enumerate() {
                acc = self.vecs[j * n + i].mul_add(uhi, acc);
            }
            *uj = acc;
        }
    }
}
