//! Analysis drivers over a sealed [`System`]: stationary equilibrium,
//! linearised harmonic response, modal analysis, periodic steady state
//! (shooting and harmonic balance via `fs-orbit`), and forward-tangent /
//! matrix-free IFT adjoint sensitivities (via `fs-adjoint`).
//!
//! Every driver operates directly on the unchanged [`System`] and its
//! two-generation checkpoint store: no second solver graph, no duplicated
//! physics equations, and no silent fallback when a component lacks a
//! required capability tier.

use core::cell::{Cell, RefCell};
use core::f64::consts::TAU;
use core::fmt;

use fs_adjoint::{AdjointReport, ift_gradient_matfree};
use fs_exec::Cx;
use fs_la::eigen_complex::{EigError, LuComplexError, eig_with_cancel, lu_complex_with_cancel};
use fs_math::c64::C64;
use fs_orbit::{
    HbAnchor, HbBudget, HbScales, OrbitError, OrbitProblem, ScaledHbOrbit, solve_hb_scaled_seeded,
};
use fs_solver::{
    FallibleNewtonError, FallibleNonlinearProblem, GmresState, LinearOp, NewtonKrylovConfig,
    solve_fallible_newton_krylov,
};

use super::{Checkpoint, CheckpointError, Phase, System, SystemError};
use crate::diagnostics::{BlockSolverPolicy, MonolithicConfig};
use crate::ledger::Sum;
use crate::solve::{factor_dense, solve_factored};
use crate::wiring::{Driver, PortRef};

/// Optional component capability tier required by an analysis driver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequiredTier {
    /// Port-map tangent (`jvp` and `vjp`) from [`fs_ports::component::PortTangent`].
    PortTangent,
    /// Monolithic stage residual (`residual`, `jvp_z`, `vjp_z`) from
    /// [`fs_ports::component::StageResidual`].
    StageResidual,
}

impl fmt::Display for RequiredTier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PortTangent => write!(f, "PortTangent"),
            Self::StageResidual => write!(f, "StageResidual"),
        }
    }
}

/// Analysis mode selectable on an unchanged [`System`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisKind {
    /// Ordinary coupled time stepping via [`System::step`].
    Transient,
    /// Held-boundary stationary equilibrium via [`System::solve_stationary`].
    Stationary,
    /// Linearised frequency-response transfer matrix via [`System::harmonic_response`].
    Harmonic,
    /// Linearised eigenmodes and damping ratios via [`System::modal_analysis`].
    Modal,
    /// Periodic steady state via Poincaré Newton shooting ([`System::solve_periodic_shooting`]).
    PeriodicShooting,
    /// Periodic steady state via harmonic balance through `fs-orbit` ([`System::solve_periodic_hb`]).
    PeriodicHarmonicBalance,
    /// Forward-tangent step sensitivity via [`System::step_sensitivity_tangent`].
    TangentSensitivity,
    /// Reverse-mode IFT adjoint sensitivity via [`System::step_sensitivity_adjoint`] or
    /// [`System::trajectory_adjoint`].
    AdjointSensitivity,
}

impl AnalysisKind {
    /// Stable lowercase identifier for diagnostics and capability refusals.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Transient => "transient",
            Self::Stationary => "stationary",
            Self::Harmonic => "harmonic",
            Self::Modal => "modal",
            Self::PeriodicShooting => "periodic-shooting",
            Self::PeriodicHarmonicBalance => "periodic-harmonic-balance",
            Self::TangentSensitivity => "tangent-sensitivity",
            Self::AdjointSensitivity => "adjoint-sensitivity",
        }
    }
}

impl fmt::Display for AnalysisKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Typed refusal from any [`System`] analysis driver.
#[derive(Debug)]
pub enum AnalysisError {
    /// A component in the system lacks a capability tier required by the
    /// requested analysis.
    MissingTier {
        /// Which analysis refused.
        analysis: AnalysisKind,
        /// Name of the first component lacking the tier.
        node: String,
        /// Required optional capability tier.
        required: RequiredTier,
    },
    /// Invalid configuration, channel specification, or dimension mismatch.
    Invalid {
        /// Human-readable reason.
        what: &'static str,
    },
    /// Iterative analysis solver did not converge within its budget.
    NotConverged {
        /// Which analysis stalled.
        analysis: AnalysisKind,
        /// Iterations spent.
        iterations: usize,
        /// Final residual norm.
        residual: f64,
    },
    /// Linear or Jacobian operator became singular or produced non-finite values.
    Singular {
        /// Which analysis hit singularity.
        analysis: AnalysisKind,
        /// Context string.
        what: &'static str,
    },
    /// Underlying coupled-step driver error.
    System(SystemError),
    /// Checkpoint capture, restore, or state write error.
    Checkpoint(CheckpointError),
    /// Periodic harmonic-balance error forwarded from `fs-orbit`.
    Orbit(OrbitError),
    /// Complex eigensolver error forwarded from `fs-la`.
    Eigen(EigError),
    /// Cancellation was requested via [`Cx`].
    Cancelled,
}

impl fmt::Display for AnalysisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTier {
                analysis,
                node,
                required,
            } => write!(
                f,
                "analysis `{analysis}` refused: component `{node}` does not offer required tier `{required}`"
            ),
            Self::Invalid { what } => write!(f, "invalid analysis input: {what}"),
            Self::NotConverged {
                analysis,
                iterations,
                residual,
            } => write!(
                f,
                "analysis `{analysis}` did not converge after {iterations} iterations (residual {residual:.3e})"
            ),
            Self::Singular { analysis, what } => {
                write!(
                    f,
                    "analysis `{analysis}` encountered singular operator: {what}"
                )
            }
            Self::System(err) => write!(f, "system step refused during analysis: {err}"),
            Self::Checkpoint(err) => write!(f, "checkpoint error during analysis: {err}"),
            Self::Orbit(err) => write!(f, "periodic orbit solver refused: {err}"),
            Self::Eigen(err) => write!(f, "complex eigensolver refused: {err:?}"),
            Self::Cancelled => write!(f, "analysis cancelled"),
        }
    }
}

impl std::error::Error for AnalysisError {}

impl From<SystemError> for AnalysisError {
    fn from(value: SystemError) -> Self {
        match value {
            SystemError::Cancelled => Self::Cancelled,
            other => Self::System(other),
        }
    }
}

impl From<CheckpointError> for AnalysisError {
    fn from(value: CheckpointError) -> Self {
        Self::Checkpoint(value)
    }
}

impl From<OrbitError> for AnalysisError {
    fn from(value: OrbitError) -> Self {
        Self::Orbit(value)
    }
}

impl From<EigError> for AnalysisError {
    fn from(value: EigError) -> Self {
        match value {
            EigError::Cancelled => Self::Cancelled,
            other => Self::Eigen(other),
        }
    }
}

/// Bidirectional map between a [`System`]'s accepted component states and a
/// flat real state vector `x in R^n` used by stationary, linearisation,
/// modal, periodic, and checkpointed trajectory-adjoint drivers.
pub trait StateMap {
    /// Dimension `n` of the flat state vector.
    fn dim(&self) -> usize;

    /// Read the system's current accepted states into `out` (length `self.dim()`).
    ///
    /// # Errors
    /// [`AnalysisError`] if a component state cannot be extracted or lengths disagree.
    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError>;

    /// Write `state` (length `self.dim()`) into the system's accepted states.
    ///
    /// # Errors
    /// [`AnalysisError`] if a state entry is out of domain or type-mismatched.
    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError>;

    /// Human-readable coordinate names (defaults to `x[0]`, `x[1]`, ...).
    fn labels(&self) -> Vec<String> {
        (0..self.dim()).map(|i| format!("x[{i}]")).collect()
    }
}

/// One scalar boundary excitation channel addressed by its [`PortRef`] and
/// component offset within that port.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryChannel {
    /// Boundary port supplied to [`crate::SystemBuilder::boundary`].
    pub port: PortRef,
    /// Zero-based scalar component within the port's input width.
    pub component: usize,
}

/// One scalar output probe channel for linearised state-space and harmonic
/// response analyses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputChannel {
    /// Continuous state coordinate `i` of the caller's [`StateMap`].
    State(usize),
    /// Step-mean scalar input of `port` at offset `component`.
    PortInput {
        /// Port to observe.
        port: PortRef,
        /// Zero-based scalar component within the port's input width.
        component: usize,
    },
    /// Step-mean scalar output of `port` at offset `component`.
    PortOutput {
        /// Port to observe.
        port: PortRef,
        /// Zero-based scalar component within the port's output width.
        component: usize,
    },
}

/// Controls for [`System::solve_stationary`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationaryConfig {
    /// Probe step size [s] used to evaluate the exact implicit-midpoint fixed
    /// point `(Phi_dt(x) - x) / dt = 0` under monolithic `StageResidual` solves.
    pub probe_dt_s: f64,
    /// Absolute Newton residual tolerance on `(Phi_dt(x) - x) / dt`.
    pub absolute_tolerance: f64,
    /// Relative Newton residual tolerance (0.0 keeps acceptance purely absolute).
    pub relative_tolerance: f64,
    /// Maximum outer Newton iterations.
    pub max_iterations: usize,
    /// Maximum Krylov restart dimension.
    pub linear_restart: usize,
    /// Maximum GMRES restart cycles per Newton step.
    pub max_linear_cycles: usize,
    /// Monolithic inner-block Newton--Krylov settings used on `StageResidual` loops.
    pub monolithic: MonolithicConfig,
    /// Relative perturbation scale for directional finite-difference Jacobian actions.
    pub fd_eps: f64,
    /// When `true`, install the converged equilibrium state into the [`System`];
    /// when `false` (default), restore the pre-solve checkpoint before returning.
    pub commit: bool,
}

impl Default for StationaryConfig {
    fn default() -> Self {
        Self {
            probe_dt_s: 1.0e-3,
            absolute_tolerance: 1.0e-10,
            relative_tolerance: 0.0,
            max_iterations: 40,
            linear_restart: 24,
            max_linear_cycles: 4,
            monolithic: MonolithicConfig::default(),
            fd_eps: 1.0e-6,
            commit: false,
        }
    }
}

/// Converged stationary operating point of a [`System`].
#[derive(Debug, Clone)]
pub struct StationarySolution {
    /// Converged equilibrium state vector `x*`.
    pub state: Vec<f64>,
    /// Equilibrium scalar port inputs at `x*`.
    pub inputs: Vec<f64>,
    /// Equilibrium scalar port outputs at `x*`.
    pub outputs: Vec<f64>,
    /// Final unscaled fixed-point residual norm `||(Phi_dt(x*) - x*) / dt||_2`.
    pub residual_norm: f64,
    /// Outer Newton iterations spent.
    pub iterations: usize,
    /// Total residual evaluations performed.
    pub residual_calls: usize,
    /// Total Jacobian-vector actions performed.
    pub jacobian_calls: usize,
}

/// Controls for [`System::linearize_state_space`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearizeConfig {
    /// Probe step size [s] for one-step implicit-midpoint linearisation.
    pub dt_s: f64,
    /// Relative perturbation step for central differences about the operating point.
    pub fd_eps: f64,
}

impl Default for LinearizeConfig {
    fn default() -> Self {
        Self {
            dt_s: 1.0e-3,
            fd_eps: 1.0e-5,
        }
    }
}

/// Discrete one-step (`M_d, B_d, C_d, D_d`) and exact bilinear/Cayley
/// continuous-time (`A, B, C, D`) state-space linearisation about an accepted
/// operating state.
///
/// All matrices are stored row-major.
#[derive(Debug, Clone)]
pub struct StateSpaceLinearization {
    /// State dimension `n`.
    pub state_dim: usize,
    /// Number of perturbed boundary input channels `m`.
    pub input_dim: usize,
    /// Number of observed output channels `p`.
    pub output_dim: usize,
    /// Effective probe step length `dt` [s] used for the discrete transition.
    pub dt_s: f64,
    /// Operating state `x_0` about which the system was linearised.
    pub operating_state: Vec<f64>,
    /// Discrete state transition Jacobian `M_d = d x_{k+1} / d x_k` (`n x n`).
    pub m_d: Vec<f64>,
    /// Discrete input Jacobian `B_d = d x_{k+1} / d u_k` (`n x m`).
    pub b_d: Vec<f64>,
    /// Discrete output state Jacobian `C_d = d y_k / d x_k` (`p x n`).
    pub c_d: Vec<f64>,
    /// Discrete output feedthrough `D_d = d y_k / d u_k` (`p x m`).
    pub d_d: Vec<f64>,
    /// Continuous-time state matrix `A = (2/dt) (I + M_d)^{-1} (M_d - I)` (`n x n`).
    pub a: Vec<f64>,
    /// Continuous-time input matrix `B = (2/dt) (I + M_d)^{-1} B_d` (`n x m`).
    pub b: Vec<f64>,
    /// Continuous-time output matrix `C` (`p x n`).
    pub c: Vec<f64>,
    /// Continuous-time feedthrough matrix `D` (`p x m`).
    pub d: Vec<f64>,
}

/// Frequency-response transfer matrices evaluated from a linearised [`System`].
#[derive(Debug, Clone)]
pub struct HarmonicResponse {
    /// Underlying continuous and discrete state-space linearisation.
    pub linearization: StateSpaceLinearization,
    /// Evaluated angular frequencies `omega` [rad/s].
    pub frequencies_rad_s: Vec<f64>,
    /// Continuous transfer matrices `H(i omega) = C (i omega I - A)^{-1} B + D`,
    /// indexed as `[freq][out * input_dim + inp]`.
    pub continuous_transfer: Vec<Vec<C64>>,
    /// Discrete step transfer matrices
    /// `H_d(e^{i omega dt}) = C_d (e^{i omega dt} I - M_d)^{-1} B_d + D_d`,
    /// indexed as `[freq][out * input_dim + inp]`.
    pub discrete_transfer: Vec<Vec<C64>>,
}

impl HarmonicResponse {
    /// Continuous transfer entry `H_{out, inp}(i omega_k)`.
    #[must_use]
    pub fn continuous(&self, freq_index: usize, out_index: usize, in_index: usize) -> C64 {
        self.continuous_transfer[freq_index][out_index * self.linearization.input_dim + in_index]
    }

    /// Discrete step transfer entry `H_{d; out, inp}(e^{i omega_k dt})`.
    #[must_use]
    pub fn discrete(&self, freq_index: usize, out_index: usize, in_index: usize) -> C64 {
        self.discrete_transfer[freq_index][out_index * self.linearization.input_dim + in_index]
    }
}

/// Controls for [`System::modal_analysis`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalConfig {
    /// Linearisation settings.
    pub linearize: LinearizeConfig,
    /// Per-window QR sweep cap forwarded to `fs-la`.
    pub max_qr_sweeps: usize,
}

impl Default for ModalConfig {
    fn default() -> Self {
        Self {
            linearize: LinearizeConfig::default(),
            max_qr_sweeps: 160,
        }
    }
}

/// One eigenmode of the linearised continuous-time state matrix `A`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalMode {
    /// Complex eigenvalue `lambda = sigma + i omega_d` [1/s].
    pub eigenvalue: C64,
    /// Damped oscillation frequency `|Im(lambda)| / (2 pi)` [Hz].
    pub frequency_hz: f64,
    /// Damped angular frequency `|Im(lambda)|` [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Undamped natural angular frequency `|lambda|` [rad/s].
    pub natural_frequency_rad_s: f64,
    /// Viscous damping ratio `-Re(lambda) / |lambda|` (0 if `|lambda| == 0`).
    pub damping_ratio: f64,
    /// Exponential decay rate `-Re(lambda)` [1/s] (positive for decaying modes).
    pub decay_rate_per_s: f64,
}

/// Modal spectrum of a linearised [`System`].
#[derive(Debug, Clone)]
pub struct ModalAnalysis {
    /// Underlying state-space linearisation.
    pub linearization: StateSpaceLinearization,
    /// Eigenmodes sorted deterministically by frequency then decay rate.
    pub modes: Vec<ModalMode>,
}

/// Controls for [`System::solve_periodic_shooting`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicShootingConfig {
    /// Period `T` [s] of the forced periodic orbit.
    pub period_s: f64,
    /// Fixed number of uniform steps per period (`dt = period_s / steps_per_period`).
    pub steps_per_period: usize,
    /// Convergence tolerance on `||Phi_T(x0) - x0||_2`.
    pub tolerance: f64,
    /// Maximum outer Newton iterations.
    pub max_iterations: usize,
    /// Relative perturbation step for the monodromy matrix `D Phi_T`.
    pub fd_eps: f64,
    /// Per-window QR sweep cap for Floquet multipliers.
    pub max_qr_sweeps: usize,
    /// When `true`, install the periodic initial state `x0*` into the [`System`].
    pub commit: bool,
}

impl Default for PeriodicShootingConfig {
    fn default() -> Self {
        Self {
            period_s: 1.0,
            steps_per_period: 100,
            tolerance: 1.0e-10,
            max_iterations: 25,
            fd_eps: 1.0e-6,
            max_qr_sweeps: 160,
            commit: false,
        }
    }
}

/// Converged periodic orbit from Poincaré Newton shooting on [`System::step`].
#[derive(Debug, Clone)]
pub struct PeriodicShootingSolution {
    /// Converged periodic initial state `x0*` satisfying `Phi_T(x0*) = x0*`.
    pub initial_state: Vec<f64>,
    /// Final closure defect `||Phi_T(x0*) - x0*||_2`.
    pub closure_residual: f64,
    /// Outer Newton iterations spent.
    pub iterations: usize,
    /// Sampled states across one full period (length `steps_per_period + 1`,
    /// starting at `x0*` and ending at `Phi_T(x0*)`).
    pub trajectory: Vec<Vec<f64>>,
    /// Row-major `n x n` monodromy matrix `D Phi_T(x0*)`.
    pub monodromy: Vec<f64>,
    /// Floquet multipliers (eigenvalues of `monodromy`), sorted deterministically.
    pub floquet_multipliers: Vec<C64>,
    /// Spectral radius `max |mu_i|` (`< 1` indicates linear orbital stability).
    pub spectral_radius: f64,
}

/// Controls for [`System::solve_periodic_hb`] through `fs-orbit`.
#[derive(Debug, Clone)]
pub struct PeriodicHbConfig {
    /// Harmonic-balance anchor (forced frequency, autonomous, or backbone).
    pub anchor: HbAnchor,
    /// Harmonic truncation and Newton budgets forwarded to `fs-orbit`.
    pub budget: HbBudget,
    /// Physical conditioning scales forwarded to `fs-orbit`.
    pub scales: HbScales,
    /// Packed initial Fourier seed (`d * (2 * N + 1)` entries, plus omega slot
    /// if autonomous/backbone). If `None`, a zero/unit seed is constructed from
    /// the current accepted state.
    pub seed: Option<Vec<f64>>,
    /// Probe step size [s] for the `O(dt^2)` Richardson vector-field extraction
    /// `f(t, x) = 2 D(dt/2) - D(dt)`.
    pub probe_dt_s: f64,
}

/// Converged harmonic-balance periodic orbit solved through `fs-orbit`.
#[derive(Debug, Clone)]
pub struct PeriodicHbSolution {
    /// Scaled harmonic-balance orbit returned by `fs-orbit`.
    pub orbit: ScaledHbOrbit,
    /// Probe step size [s] used for Richardson vector-field evaluation.
    pub probe_dt_s: f64,
}

/// Controls for matrix-free tangent and IFT adjoint sensitivity solves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SensitivityConfig {
    /// Relative GMRES tolerance passed to `fs-solver` / `fs-adjoint::ift_gradient_matfree`.
    pub gmres_tolerance: f64,
    /// Maximum GMRES restart cycles.
    pub max_gmres_cycles: usize,
    /// Relative perturbation scale used only for inter-step state-transition
    /// Jacobians in multi-step [`System::trajectory_adjoint`].
    pub state_fd_eps: f64,
    /// When `true` on single-step sensitivity drivers, keep the accepted step
    /// published; when `false` (default), restore the pre-step checkpoint.
    pub commit_step: bool,
}

impl Default for SensitivityConfig {
    fn default() -> Self {
        Self {
            gmres_tolerance: 1.0e-12,
            max_gmres_cycles: 20,
            state_fd_eps: 1.0e-5,
            commit_step: false,
        }
    }
}

/// Partial derivatives of a smooth scalar objective `J(u, y, b)` on one step.
#[derive(Debug, Clone)]
pub struct StepObjective {
    /// Partial derivative `dJ / d inputs` (length `system.input_len()`, or empty for zero).
    pub dj_du: Vec<f64>,
    /// Partial derivative `dJ / d outputs` (length `system.output_len()`, or empty for zero).
    pub dj_dy: Vec<f64>,
    /// Explicit partial derivative `dJ / d boundaries` (length `system.boundary_scalar_count()`,
    /// or empty for zero).
    pub dj_db: Vec<f64>,
}

/// Forward-tangent step sensitivity along a boundary perturbation `db`.
#[derive(Debug, Clone)]
pub struct StepTangentSensitivity {
    /// Effective step length taken [s].
    pub dt_s: f64,
    /// Linearised change in all scalar port inputs `du`.
    pub d_inputs: Vec<f64>,
    /// Linearised change in all scalar port outputs `dy`.
    pub d_outputs: Vec<f64>,
    /// Forward GMRES iteration count.
    pub gmres_iters: usize,
    /// Achieved GMRES relative residual.
    pub gmres_residual: f64,
}

/// Reverse-mode IFT adjoint sensitivity of one coupled step.
#[derive(Debug, Clone)]
pub struct StepAdjointSensitivity {
    /// Effective step length taken [s].
    pub dt_s: f64,
    /// Total gradient `dJ / d b` for every scalar boundary channel in declaration order.
    pub boundary_gradient: Vec<f64>,
    /// Optional custom parameter gradient `dJ / d p` when a parameter pullback
    /// was supplied to [`System::step_sensitivity_adjoint_with_pullback`].
    pub parameter_gradient: Vec<f64>,
    /// Interface adjoint multiplier `lambda` solving `(dR/du)^T lambda = dJ_total/du`.
    pub adjoint_inputs: Vec<f64>,
    /// Total output cotangent `bar_y = dJ/dy + (dD/dy)^T lambda` suitable for
    /// pulling back into card-derived component parameters via `(dy/dp)^T bar_y`.
    pub adjoint_outputs: Vec<f64>,
    /// Solver telemetry returned directly by `fs_adjoint::ift_gradient_matfree`.
    pub report: AdjointReport,
}

/// Per-step objective partials for checkpointed multi-step trajectory adjoints.
#[derive(Debug, Clone, Default)]
pub struct TrajectoryStepPartials {
    /// Partial derivative `dJ / d state_k` at the start of step `k` (empty for zero).
    pub dj_ds: Vec<f64>,
    /// Partial derivative `dJ / d inputs_k` (empty for zero).
    pub dj_du: Vec<f64>,
    /// Partial derivative `dJ / d outputs_k` (empty for zero).
    pub dj_dy: Vec<f64>,
    /// Explicit partial derivative `dJ / d boundaries_k` (empty for zero).
    pub dj_db: Vec<f64>,
}

/// Result of a checkpointed reverse-mode trajectory adjoint sweep.
#[derive(Debug, Clone)]
pub struct TrajectoryAdjointResult {
    /// Number of steps reversed.
    pub steps: usize,
    /// Exact step lengths taken along the primal trajectory [s].
    pub step_dts_s: Vec<f64>,
    /// Terminal state `x_N` reached at the end of the primal rollout.
    pub final_state: Vec<f64>,
    /// Total gradient `dJ / d x_0` with respect to the initial accepted state.
    pub grad_initial_state: Vec<f64>,
    /// Total gradient `dJ / d b_k` with respect to each step's boundary vector.
    pub grad_boundary_per_step: Vec<Vec<f64>>,
    /// Total gradient `dJ / d b` for a constant boundary offset applied across all steps.
    pub grad_boundary_total: Vec<f64>,
    /// Per-step IFT adjoint solver reports from `fs_adjoint::ift_gradient_matfree`.
    pub step_reports: Vec<AdjointReport>,
}

/// Unified analysis dispatch request so one entry point switches analysis on an
/// unchanged [`System`].
pub enum AnalysisRequest<'a> {
    /// Run `steps` transient steps of size `dt_s`.
    Transient {
        /// Step length [s].
        dt_s: f64,
        /// Number of steps to take.
        steps: usize,
    },
    /// Solve stationary equilibrium via `StageResidual`.
    Stationary {
        /// State coordinate map.
        state_map: &'a dyn StateMap,
        /// Stationary solver configuration.
        config: StationaryConfig,
    },
    /// Evaluate linearised harmonic response transfer matrices.
    Harmonic {
        /// State coordinate map.
        state_map: &'a dyn StateMap,
        /// Boundary excitation channels.
        inputs: &'a [BoundaryChannel],
        /// Output probe channels.
        outputs: &'a [OutputChannel],
        /// Angular frequencies `omega` [rad/s].
        frequencies_rad_s: &'a [f64],
        /// Linearisation controls.
        config: LinearizeConfig,
    },
    /// Compute linearised modal spectrum.
    Modal {
        /// State coordinate map.
        state_map: &'a dyn StateMap,
        /// Modal solver controls.
        config: ModalConfig,
    },
    /// Solve periodic steady state via Poincaré Newton shooting.
    PeriodicShooting {
        /// State coordinate map.
        state_map: &'a dyn StateMap,
        /// Optional initial state guess.
        initial_guess: Option<&'a [f64]>,
        /// Shooting controls.
        config: PeriodicShootingConfig,
    },
    /// Solve periodic steady state via `fs-orbit` harmonic balance.
    PeriodicHarmonicBalance {
        /// State coordinate map.
        state_map: &'a dyn StateMap,
        /// Harmonic-balance controls.
        config: PeriodicHbConfig,
    },
    /// Compute single-step matrix-free forward tangent sensitivity.
    TangentSensitivity {
        /// Step length [s].
        dt_s: f64,
        /// Boundary perturbation vector (length [`System::boundary_scalar_count`]).
        d_boundaries: &'a [f64],
        /// Sensitivity solver controls.
        config: SensitivityConfig,
    },
    /// Compute single-step matrix-free IFT adjoint sensitivity.
    AdjointSensitivity {
        /// Step length [s].
        dt_s: f64,
        /// Step objective partial derivatives.
        objective: StepObjective,
        /// Sensitivity solver controls.
        config: SensitivityConfig,
    },
}

/// Typed outcome returned by [`System::run_analysis`].
#[derive(Debug)]
pub enum AnalysisOutcome {
    /// Completed transient step reports.
    Transient(Vec<super::StepReport>),
    /// Converged stationary operating point.
    Stationary(StationarySolution),
    /// Linearised harmonic frequency response.
    Harmonic(HarmonicResponse),
    /// Linearised modal spectrum.
    Modal(ModalAnalysis),
    /// Converged shooting periodic orbit.
    PeriodicShooting(PeriodicShootingSolution),
    /// Converged harmonic-balance periodic orbit from `fs-orbit`.
    PeriodicHarmonicBalance(PeriodicHbSolution),
    /// Single-step forward tangent sensitivity.
    TangentSensitivity(StepTangentSensitivity),
    /// Single-step IFT adjoint sensitivity.
    AdjointSensitivity(StepAdjointSensitivity),
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

struct SystemInterfaceOp<'a, 'cx> {
    system: &'a System,
    old_states: &'a [Box<dyn core::any::Any>],
    dt_s: f64,
    inputs: &'a [f64],
    outputs: &'a [f64],
    cx: &'a Cx<'cx>,
    cancelled: Cell<bool>,
    error: RefCell<Option<AnalysisError>>,
}

impl<'a, 'cx> SystemInterfaceOp<'a, 'cx> {
    fn new(
        system: &'a System,
        old_states: &'a [Box<dyn core::any::Any>],
        dt_s: f64,
        inputs: &'a [f64],
        outputs: &'a [f64],
        cx: &'a Cx<'cx>,
    ) -> Self {
        Self {
            system,
            old_states,
            dt_s,
            inputs,
            outputs,
            cx,
            cancelled: Cell::new(false),
            error: RefCell::new(None),
        }
    }

    fn latch_refusal(&self, node: usize, error: crate::erased::OwnerError) {
        if self.cx.is_cancel_requested() {
            self.cancelled.set(true);
        } else if self.error.borrow().is_none() {
            *self.error.borrow_mut() = Some(AnalysisError::System(SystemError::Component {
                node: self.system.names[node].clone(),
                phase: Phase::Propose,
                error,
            }));
        }
    }

    fn check_status(&self) -> Result<(), AnalysisError> {
        if self.cancelled.get() || self.cx.is_cancel_requested() {
            return Err(AnalysisError::Cancelled);
        }
        if let Some(err) = self.error.borrow_mut().take() {
            return Err(err);
        }
        Ok(())
    }

    /// Evaluate `dy = (df / du) v` across all nodes using exact `PortTangent::jvp`.
    fn forward_port_jvp(&self, v: &[f64], dy: &mut [f64]) -> Result<(), AnalysisError> {
        dy.fill(0.0);
        for node in 0..self.system.nodes.len() {
            if self.cx.is_cancel_requested() {
                self.cancelled.set(true);
                return Err(AnalysisError::Cancelled);
            }
            let inp_r = self.system.node_inputs(node);
            let out_r = self.system.node_outputs(node);
            if inp_r.is_empty() || out_r.is_empty() {
                continue;
            }
            self.system.nodes[node]
                .jvp(
                    self.old_states[node].as_ref(),
                    self.dt_s,
                    &self.inputs[inp_r.clone()],
                    &v[inp_r],
                    &mut dy[out_r],
                    self.cx,
                )
                .map_err(|err| {
                    if self.cx.is_cancel_requested() {
                        AnalysisError::Cancelled
                    } else {
                        AnalysisError::System(SystemError::Component {
                            node: self.system.names[node].clone(),
                            phase: Phase::Propose,
                            error: err,
                        })
                    }
                })?;
        }
        Ok(())
    }

    /// Evaluate `bar_u = (df / du)^T bar_y` across all nodes using exact `PortTangent::vjp`.
    fn reverse_port_vjp(&self, bar_y: &[f64], bar_u: &mut [f64]) -> Result<(), AnalysisError> {
        bar_u.fill(0.0);
        for node in 0..self.system.nodes.len() {
            if self.cx.is_cancel_requested() {
                self.cancelled.set(true);
                return Err(AnalysisError::Cancelled);
            }
            let inp_r = self.system.node_inputs(node);
            let out_r = self.system.node_outputs(node);
            if inp_r.is_empty() || out_r.is_empty() {
                continue;
            }
            self.system.nodes[node]
                .vjp(
                    self.old_states[node].as_ref(),
                    self.dt_s,
                    &self.inputs[inp_r.clone()],
                    &bar_y[out_r],
                    &mut bar_u[inp_r],
                    self.cx,
                )
                .map_err(|err| {
                    if self.cx.is_cancel_requested() {
                        AnalysisError::Cancelled
                    } else {
                        AnalysisError::System(SystemError::Component {
                            node: self.system.names[node].clone(),
                            phase: Phase::Propose,
                            error: err,
                        })
                    }
                })?;
        }
        Ok(())
    }

    /// Scatter `bar_y += (dD / dy)^T w`.
    fn scatter_driver_transpose(&self, w: &[f64], bar_y: &mut [f64]) {
        let mut grad = Vec::new();
        for (i, &wi) in w.iter().enumerate() {
            if wi == 0.0 {
                continue;
            }
            self.system.drivers[i].gradient(self.outputs, &mut grad);
            for &(src, slope) in &grad {
                bar_y[src] += slope * wi;
            }
        }
    }
}

impl LinearOp for SystemInterfaceOp<'_, '_> {
    fn n(&self) -> usize {
        self.inputs.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.cancelled.get() || self.cx.is_cancel_requested() {
            self.cancelled.set(true);
            y.fill(f64::NAN);
            return;
        }
        if self.error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        let mut dy = vec![0.0_f64; self.outputs.len()];
        for node in 0..self.system.nodes.len() {
            let inp_r = self.system.node_inputs(node);
            let out_r = self.system.node_outputs(node);
            if inp_r.is_empty() || out_r.is_empty() {
                continue;
            }
            if let Err(err) = self.system.nodes[node].jvp(
                self.old_states[node].as_ref(),
                self.dt_s,
                &self.inputs[inp_r.clone()],
                &x[inp_r],
                &mut dy[out_r],
                self.cx,
            ) {
                self.latch_refusal(node, err);
                y.fill(f64::NAN);
                return;
            }
        }
        let mut grad = Vec::new();
        for i in 0..self.inputs.len() {
            self.system.drivers[i].gradient(self.outputs, &mut grad);
            let mut d_target = 0.0_f64;
            for &(src, slope) in &grad {
                d_target += slope * dy[src];
            }
            y[i] = x[i] - d_target;
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        if self.cancelled.get() || self.cx.is_cancel_requested() {
            self.cancelled.set(true);
            y.fill(f64::NAN);
            return;
        }
        if self.error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        let mut bar_y = vec![0.0_f64; self.outputs.len()];
        self.scatter_driver_transpose(x, &mut bar_y);
        let mut bar_u = vec![0.0_f64; self.inputs.len()];
        for node in 0..self.system.nodes.len() {
            let inp_r = self.system.node_inputs(node);
            let out_r = self.system.node_outputs(node);
            if inp_r.is_empty() || out_r.is_empty() {
                continue;
            }
            if let Err(err) = self.system.nodes[node].vjp(
                self.old_states[node].as_ref(),
                self.dt_s,
                &self.inputs[inp_r.clone()],
                &bar_y[out_r],
                &mut bar_u[inp_r],
                self.cx,
            ) {
                self.latch_refusal(node, err);
                y.fill(f64::NAN);
                return;
            }
        }
        for i in 0..self.inputs.len() {
            y[i] = x[i] - bar_u[i];
        }
    }
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn clamp_perturbation(base: f64, eps: f64) -> f64 {
    let raw = eps * base.abs().max(1.0);
    if base > 0.0 && raw >= 0.25 * base {
        0.25 * base
    } else {
        raw
    }
}

type LinearizedMatrices = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

/// Compute continuous-time `(A, B, C, D)` from discrete implicit-midpoint
/// `(M_d, B_d, C_d, D_d)` via exact Cayley/bilinear inversion.
#[allow(clippy::too_many_arguments)]
fn cayley_to_continuous(
    analysis: AnalysisKind,
    n: usize,
    m: usize,
    p: usize,
    dt_s: f64,
    output_channels: &[OutputChannel],
    m_d: &[f64],
    b_d: &[f64],
    c_d: &[f64],
    d_d: &[f64],
) -> Result<LinearizedMatrices, AnalysisError> {
    // Factor P_mat = I + M_d once (n x n).
    let mut p_factors = vec![0.0_f64; n * n];
    for r in 0..n {
        for c in 0..n {
            p_factors[r * n + c] = m_d[r * n + c] + if r == c { 1.0 } else { 0.0 };
        }
    }
    let mut pivots = vec![0_usize; n];
    if !factor_dense(&mut p_factors, &mut pivots, n) {
        return Err(AnalysisError::Singular {
            analysis,
            what: "I + M_d is singular (Nyquist pole at z = -1); try a smaller dt_s",
        });
    }

    // A = (2 / dt) (I + M_d)^{-1} (M_d - I), column by column.
    let scale = 2.0 / dt_s;
    let mut a = vec![0.0_f64; n * n];
    let mut col = vec![0.0_f64; n];
    for j in 0..n {
        for i in 0..n {
            col[i] = scale * (m_d[i * n + j] - if i == j { 1.0 } else { 0.0 });
        }
        if !solve_factored(&p_factors, &pivots, &mut col, n) {
            return Err(AnalysisError::Singular {
                analysis,
                what: "Cayley solve for A produced non-finite values",
            });
        }
        for i in 0..n {
            a[i * n + j] = col[i];
        }
    }

    // B = (2 / dt) (I + M_d)^{-1} B_d.
    let mut b = vec![0.0_f64; n * m];
    for j in 0..m {
        for i in 0..n {
            col[i] = scale * b_d[i * m + j];
        }
        if !solve_factored(&p_factors, &pivots, &mut col, n) {
            return Err(AnalysisError::Singular {
                analysis,
                what: "Cayley solve for B produced non-finite values",
            });
        }
        for i in 0..n {
            b[i * m + j] = col[i];
        }
    }

    // Factor (I + M_d)^T to solve C_row = 2 C_d_row (I + M_d)^{-1} for step-mean ports.
    let mut pt_factors = vec![0.0_f64; n * n];
    for r in 0..n {
        for c in 0..n {
            pt_factors[r * n + c] = m_d[c * n + r] + if r == c { 1.0 } else { 0.0 };
        }
    }
    let mut pt_pivots = vec![0_usize; n];
    if !factor_dense(&mut pt_factors, &mut pt_pivots, n) {
        return Err(AnalysisError::Singular {
            analysis,
            what: "(I + M_d)^T is singular",
        });
    }

    let mut c = vec![0.0_f64; p * n];
    let mut d = vec![0.0_f64; p * m];
    for (row, &channel) in output_channels.iter().enumerate() {
        match channel {
            OutputChannel::State(idx) => {
                c[row * n + idx] = 1.0;
                // Continuous feedthrough for a pure state coordinate is identically zero.
            }
            OutputChannel::PortInput { .. } | OutputChannel::PortOutput { .. } => {
                for j in 0..n {
                    col[j] = 2.0 * c_d[row * n + j];
                }
                if !solve_factored(&pt_factors, &pt_pivots, &mut col, n) {
                    return Err(AnalysisError::Singular {
                        analysis,
                        what: "Cayley solve for C produced non-finite values",
                    });
                }
                for j in 0..n {
                    c[row * n + j] = col[j];
                }
                for k in 0..m {
                    let mut c_bd = 0.0_f64;
                    for j in 0..n {
                        c_bd += c[row * n + j] * b_d[j * m + k];
                    }
                    d[row * m + k] = d_d[row * m + k] - 0.5 * c_bd;
                }
            }
        }
    }
    Ok((a, b, c, d))
}

impl System {
    // -----------------------------------------------------------------------
    // Indexing and controlled state/boundary accessors for analysis drivers
    // -----------------------------------------------------------------------

    /// Total number of scalar port inputs across all components in the system.
    #[must_use]
    pub fn input_len(&self) -> usize {
        self.inputs.len()
    }

    /// Total number of scalar port outputs across all components in the system.
    #[must_use]
    pub fn output_len(&self) -> usize {
        self.outputs.len()
    }

    /// Total number of scalar boundary input entries in declaration order.
    #[must_use]
    pub fn boundary_scalar_count(&self) -> usize {
        self.boundary_values.iter().map(Vec::len).sum()
    }

    /// All scalar boundary channels in declaration order.
    #[must_use]
    pub fn boundary_channels(&self) -> Vec<BoundaryChannel> {
        let mut channels = Vec::with_capacity(self.boundary_scalar_count());
        for (port, vals) in self
            .boundaries
            .iter()
            .map(|(p, _)| *p)
            .zip(&self.boundary_values)
        {
            for component in 0..vals.len() {
                channels.push(BoundaryChannel { port, component });
            }
        }
        channels
    }

    /// Flat scalar index into [`Self::input_len`] for `(port, component)`.
    #[must_use]
    pub fn input_index(&self, port: PortRef, component: usize) -> Option<usize> {
        let node = self.layouts.get(port.node)?;
        if port.port >= self.nodes[port.node].ports().len() {
            return None;
        }
        let range = node.input_range(port.port);
        (component < range.len()).then_some(self.input_start[port.node] + range.start + component)
    }

    /// Flat scalar index into [`Self::output_len`] for `(port, component)`.
    #[must_use]
    pub fn output_index(&self, port: PortRef, component: usize) -> Option<usize> {
        let node = self.layouts.get(port.node)?;
        if port.port >= self.nodes[port.node].ports().len() {
            return None;
        }
        let range = node.output_range(port.port);
        (component < range.len()).then_some(self.output_start[port.node] + range.start + component)
    }

    /// Replace the accepted state of `node`, advancing the accepted generation
    /// and invalidating internal Jacobian caches.
    ///
    /// # Errors
    /// [`CheckpointError::State`] if `node` is out of range or `S` does not
    /// match the component's concrete state type.
    pub fn set_accepted_state<S: 'static>(
        &mut self,
        node: crate::wiring::NodeId,
        value: S,
    ) -> Result<(), CheckpointError> {
        let Some(slot) = self.accepted.get_mut(node.0) else {
            return Err(CheckpointError::Image);
        };
        let Some(typed) = slot.downcast_mut::<S>() else {
            return Err(CheckpointError::State {
                node: self.names[node.0].clone(),
                source: self.nodes[node.0]
                    .clone_state(&())
                    .err()
                    .unwrap_or_else(|| {
                        // Fallback unreachable unless state type is `()` which would have matched.
                        self.nodes[node.0]
                            .clone_state(&0_u8)
                            .expect_err("type mismatch")
                    }),
            });
        };
        *typed = value;
        self.accepted_generation = self.accepted_generation.wrapping_add(1);
        self.invalidate_caches();
        Ok(())
    }

    /// Set a constant additive bias on boundary channel `(port, component)`
    /// (applied after evaluating its [`super::BoundaryFn`] during trials).
    /// Useful for perturbations and independent gradient checks on the exact
    /// same sealed [`System`].
    ///
    /// # Errors
    /// [`AnalysisError::Invalid`] if `(port, component)` is not a declared boundary
    /// or `bias` is non-finite.
    pub fn set_boundary_bias(
        &mut self,
        port: PortRef,
        component: usize,
        bias: f64,
    ) -> Result<(), AnalysisError> {
        if !bias.is_finite() {
            return Err(AnalysisError::Invalid {
                what: "boundary bias must be finite",
            });
        }
        let (b_idx, comp) = self.resolve_boundary_channel(BoundaryChannel { port, component })?;
        self.boundary_bias[b_idx][comp] = bias;
        Ok(())
    }

    /// Reset all additive boundary biases to zero.
    pub fn clear_boundary_biases(&mut self) {
        for slot in &mut self.boundary_bias {
            slot.fill(0.0);
        }
    }

    fn resolve_boundary_channel(
        &self,
        channel: BoundaryChannel,
    ) -> Result<(usize, usize), AnalysisError> {
        let Some(b_idx) = self
            .boundaries
            .iter()
            .position(|(port, _)| *port == channel.port)
        else {
            return Err(AnalysisError::Invalid {
                what: "boundary channel port is not a declared system boundary",
            });
        };
        if channel.component >= self.boundary_values[b_idx].len() {
            return Err(AnalysisError::Invalid {
                what: "boundary channel component index out of range",
            });
        }
        Ok((b_idx, channel.component))
    }

    fn flat_boundary_map(&self) -> Vec<(usize, usize, usize)> {
        // Each entry is `(boundary_idx, component_idx, global_input_idx)`.
        let mut map = Vec::with_capacity(self.boundary_scalar_count());
        for (b_idx, vals) in self.boundary_values.iter().enumerate() {
            for comp in 0..vals.len() {
                let input_idx = self
                    .drivers
                    .iter()
                    .position(|d| {
                        matches!(
                            d,
                            Driver::Boundary {
                                boundary,
                                component
                            } if *boundary == b_idx && *component == comp
                        )
                    })
                    .expect("every declared boundary scalar drives one system input");
                map.push((b_idx, comp, input_idx));
            }
        }
        map
    }

    // -----------------------------------------------------------------------
    // Capability gating & unified switch
    // -----------------------------------------------------------------------

    /// Verify that every component in this [`System`] offers the capability
    /// tiers required by `kind`, refusing with the exact component name and
    /// missing tier otherwise.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`].
    pub fn check_analysis_capabilities(&self, kind: AnalysisKind) -> Result<(), AnalysisError> {
        match kind {
            AnalysisKind::Transient => Ok(()),
            AnalysisKind::Stationary => {
                for (idx, node) in self.nodes.iter().enumerate() {
                    if !node.has_stage_residual() {
                        return Err(AnalysisError::MissingTier {
                            analysis: kind,
                            node: self.names[idx].clone(),
                            required: RequiredTier::StageResidual,
                        });
                    }
                }
                Ok(())
            }
            AnalysisKind::Harmonic
            | AnalysisKind::Modal
            | AnalysisKind::PeriodicShooting
            | AnalysisKind::PeriodicHarmonicBalance
            | AnalysisKind::TangentSensitivity
            | AnalysisKind::AdjointSensitivity => {
                for (idx, node) in self.nodes.iter().enumerate() {
                    if !node.has_tangent() {
                        return Err(AnalysisError::MissingTier {
                            analysis: kind,
                            node: self.names[idx].clone(),
                            required: RequiredTier::PortTangent,
                        });
                    }
                }
                Ok(())
            }
        }
    }

    /// Run any supported analysis on this unchanged [`System`] from one unified
    /// request enum.
    ///
    /// # Errors
    /// Any [`AnalysisError`], including capability-gating refusals.
    pub fn run_analysis(
        &mut self,
        request: AnalysisRequest<'_>,
        cx: &Cx<'_>,
    ) -> Result<AnalysisOutcome, AnalysisError> {
        match request {
            AnalysisRequest::Transient { dt_s, steps } => {
                let mut reports = Vec::with_capacity(steps);
                for _ in 0..steps {
                    reports.push(self.step(dt_s, cx)?);
                }
                Ok(AnalysisOutcome::Transient(reports))
            }
            AnalysisRequest::Stationary { state_map, config } => Ok(AnalysisOutcome::Stationary(
                self.solve_stationary(state_map, config, cx)?,
            )),
            AnalysisRequest::Harmonic {
                state_map,
                inputs,
                outputs,
                frequencies_rad_s,
                config,
            } => Ok(AnalysisOutcome::Harmonic(self.harmonic_response(
                state_map,
                inputs,
                outputs,
                frequencies_rad_s,
                config,
                cx,
            )?)),
            AnalysisRequest::Modal { state_map, config } => Ok(AnalysisOutcome::Modal(
                self.modal_analysis(state_map, config, cx)?,
            )),
            AnalysisRequest::PeriodicShooting {
                state_map,
                initial_guess,
                config,
            } => Ok(AnalysisOutcome::PeriodicShooting(
                self.solve_periodic_shooting(state_map, initial_guess, config, cx)?,
            )),
            AnalysisRequest::PeriodicHarmonicBalance { state_map, config } => {
                Ok(AnalysisOutcome::PeriodicHarmonicBalance(
                    self.solve_periodic_hb(state_map, config, cx)?,
                ))
            }
            AnalysisRequest::TangentSensitivity {
                dt_s,
                d_boundaries,
                config,
            } => Ok(AnalysisOutcome::TangentSensitivity(
                self.step_sensitivity_tangent(dt_s, d_boundaries, config, cx)?,
            )),
            AnalysisRequest::AdjointSensitivity {
                dt_s,
                objective,
                config,
            } => Ok(AnalysisOutcome::AdjointSensitivity(
                self.step_sensitivity_adjoint(dt_s, &objective, config, cx)?,
            )),
        }
    }

    // -----------------------------------------------------------------------
    // 1. Stationary equilibrium over StageResidual blocks
    // -----------------------------------------------------------------------

    /// Solve for a held-boundary stationary equilibrium `Phi_dt(x*) = x*` using
    /// `fs-solver`'s fallible Newton--Krylov driver and monolithic
    /// [`fs_ports::component::StageResidual`] block solves.
    ///
    /// Because `Phi_dt(x) = x` holds identically if and only if the continuous
    /// balance vanishes for every `dt > 0`, this solves both damped and
    /// lossless systems directly without relaxation-to-rest marching.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `StageResidual`;
    /// [`AnalysisError::Invalid`], [`AnalysisError::NotConverged`],
    /// [`AnalysisError::System`], or [`AnalysisError::Cancelled`]. On error or
    /// when `!config.commit`, the system is restored to its pre-call checkpoint.
    #[allow(clippy::too_many_lines, clippy::items_after_statements)]
    pub fn solve_stationary<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        config: StationaryConfig,
        cx: &Cx<'_>,
    ) -> Result<StationarySolution, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::Stationary)?;
        let n = state_map.dim();
        if n == 0
            || !config.probe_dt_s.is_finite()
            || config.probe_dt_s <= 0.0
            || !config.absolute_tolerance.is_finite()
            || config.absolute_tolerance <= 0.0
            || !config.relative_tolerance.is_finite()
            || config.relative_tolerance < 0.0
            || config.max_iterations == 0
            || !config.fd_eps.is_finite()
            || config.fd_eps <= 0.0
        {
            return Err(AnalysisError::Invalid {
                what: "stationary solve requires positive state dim, probe_dt_s, tolerances, fd_eps, and max_iterations",
            });
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let checkpoint = self.checkpoint()?;
        let saved_policy = self.solver_policy;
        let mut x0 = vec![0.0_f64; n];
        state_map.read(self, &mut x0)?;

        // Use monolithic StageResidual solve on looped blocks whenever every
        // looped node declares at least one stage unknown.
        let all_looped_have_unknowns = self.blocks.iter().filter(|b| b.looped()).all(|b| {
            b.nodes
                .iter()
                .all(|&node| self.nodes[node].stage_unknowns() > 0)
        });
        if all_looped_have_unknowns {
            self.set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(config.monolithic));
        }

        struct StationaryState<'a, S: StateMap + ?Sized> {
            system: RefCell<&'a mut System>,
            checkpoint: &'a Checkpoint,
            state_map: &'a S,
            dim: usize,
            probe_dt_s: f64,
            fd_eps: f64,
            residual_calls: Cell<usize>,
            jacobian_calls: Cell<usize>,
        }

        impl<S: StateMap + ?Sized> StationaryState<'_, S> {
            fn eval_fixed_point_rate(
                &self,
                x: &[f64],
                out: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), AnalysisError> {
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                if !x.iter().all(|v| v.is_finite()) {
                    return Err(AnalysisError::Singular {
                        analysis: AnalysisKind::Stationary,
                        what: "non-finite state iterate",
                    });
                }
                let mut sys = self.system.borrow_mut();
                sys.restore(self.checkpoint)?;
                self.state_map.write(&mut sys, x)?;
                let report = sys.step(self.probe_dt_s, cx)?;
                let mut x_next = vec![0.0_f64; self.dim];
                self.state_map.read(&sys, &mut x_next)?;
                for i in 0..self.dim {
                    out[i] = (x_next[i] - x[i]) / report.dt_s;
                }
                Ok(())
            }
        }

        impl<S: StateMap + ?Sized> FallibleNonlinearProblem for StationaryState<'_, S> {
            type Error = AnalysisError;

            fn dimension(&self) -> usize {
                self.dim
            }

            fn residual(
                &self,
                x: &[f64],
                residual: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), AnalysisError> {
                self.residual_calls.set(self.residual_calls.get() + 1);
                self.eval_fixed_point_rate(x, residual, cx)
            }

            fn jacobian_apply(
                &self,
                x: &[f64],
                direction: &[f64],
                output: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), AnalysisError> {
                self.jacobian_calls.set(self.jacobian_calls.get() + 1);
                let d_norm = l2_norm(direction);
                if d_norm == 0.0 {
                    output.fill(0.0);
                    return Ok(());
                }
                let x_scale = x.iter().map(|v| v.abs()).fold(1.0_f64, f64::max);
                let eps = self.fd_eps * x_scale / d_norm;
                let mut xp = vec![0.0_f64; self.dim];
                let mut xm = vec![0.0_f64; self.dim];
                for i in 0..self.dim {
                    xp[i] = x[i] + eps * direction[i];
                    xm[i] = x[i] - eps * direction[i];
                }
                let mut fp = vec![0.0_f64; self.dim];
                let mut fm = vec![0.0_f64; self.dim];
                match (
                    self.eval_fixed_point_rate(&xp, &mut fp, cx),
                    self.eval_fixed_point_rate(&xm, &mut fm, cx),
                ) {
                    (Ok(()), Ok(())) => {
                        let inv_2eps = 0.5 / eps;
                        for i in 0..self.dim {
                            output[i] = (fp[i] - fm[i]) * inv_2eps;
                        }
                        Ok(())
                    }
                    (Err(AnalysisError::Cancelled), _) | (_, Err(AnalysisError::Cancelled)) => {
                        Err(AnalysisError::Cancelled)
                    }
                    (Ok(()), Err(_)) => {
                        // One-sided fallback if `-eps` stepped outside a positive state domain.
                        let mut f0 = vec![0.0_f64; self.dim];
                        self.eval_fixed_point_rate(x, &mut f0, cx)?;
                        for i in 0..self.dim {
                            output[i] = (fp[i] - f0[i]) / eps;
                        }
                        Ok(())
                    }
                    (Err(_), Ok(())) => {
                        let mut f0 = vec![0.0_f64; self.dim];
                        self.eval_fixed_point_rate(x, &mut f0, cx)?;
                        for i in 0..self.dim {
                            output[i] = (f0[i] - fm[i]) / eps;
                        }
                        Ok(())
                    }
                    (Err(e), Err(_)) => Err(e),
                }
            }
        }

        let nk_cfg = NewtonKrylovConfig {
            absolute_tolerance: config.absolute_tolerance,
            relative_tolerance: config.relative_tolerance,
            linear_restart: config.linear_restart.max(n).min(64),
            max_linear_cycles: config.max_linear_cycles,
            // Central-FD JFNK on `(Phi_dt(x) - x) / dt` has linearity floor
            // `O(eps_mach / (dt * fd_eps))`; keep the Eisenstat-Walker forcing
            // floor above that cancellation noise while still achieving >= 1e4
            // residual reduction per Newton step.
            forcing_minimum: 1.0e-4,
            forcing_maximum: 1.0e-3,
            ..NewtonKrylovConfig::default()
        };

        let (solved, residual_calls, jacobian_calls) = {
            let problem = StationaryState {
                system: RefCell::new(self),
                checkpoint: &checkpoint,
                state_map,
                dim: n,
                probe_dt_s: config.probe_dt_s,
                fd_eps: config.fd_eps,
                residual_calls: Cell::new(0),
                jacobian_calls: Cell::new(0),
            };
            let solved =
                solve_fallible_newton_krylov(&problem, x0, nk_cfg, config.max_iterations, cx);
            (
                solved,
                problem.residual_calls.get(),
                problem.jacobian_calls.get(),
            )
        };

        // Always restore policy and checkpoint first; commit explicitly only on success.
        self.set_solver_policy(saved_policy);
        let sol = match solved {
            Ok(sol) => sol,
            Err(FallibleNewtonError::Cancelled) => {
                let _ = self.restore(&checkpoint);
                return Err(AnalysisError::Cancelled);
            }
            Err(FallibleNewtonError::Problem(err)) => {
                let _ = self.restore(&checkpoint);
                return Err(err);
            }
            Err(_) => {
                let _ = self.restore(&checkpoint);
                return Err(AnalysisError::Singular {
                    analysis: AnalysisKind::Stationary,
                    what: "Newton-Krylov linear solve failed",
                });
            }
        };

        let res_norm = l2_norm(&sol.residual);
        if !sol.report.converged || !res_norm.is_finite() {
            let _ = self.restore(&checkpoint);
            return Err(AnalysisError::NotConverged {
                analysis: AnalysisKind::Stationary,
                iterations: sol.report.iterations,
                residual: res_norm,
            });
        }

        // Capture equilibrium port values at the converged operating state.
        self.restore(&checkpoint)?;
        state_map.write(self, &sol.x)?;
        if all_looped_have_unknowns {
            self.set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(config.monolithic));
        }
        let step_res = self.step(config.probe_dt_s, cx);
        self.set_solver_policy(saved_policy);
        step_res?;
        let eq_inputs = self.accepted_inputs.clone();
        let eq_outputs = self.accepted_outputs.clone();

        // Restore pre-solve checkpoint; optionally install equilibrium state at t0.
        self.restore(&checkpoint)?;
        if config.commit {
            state_map.write(self, &sol.x)?;
            self.accepted_inputs.copy_from_slice(&eq_inputs);
            self.accepted_outputs.copy_from_slice(&eq_outputs);
        }

        Ok(StationarySolution {
            state: sol.x,
            inputs: eq_inputs,
            outputs: eq_outputs,
            residual_norm: res_norm,
            iterations: sol.report.iterations,
            residual_calls,
            jacobian_calls,
        })
    }

    // -----------------------------------------------------------------------
    // 2. Linearised state-space, harmonic response, and modal analysis
    // -----------------------------------------------------------------------

    fn read_output_channels<S: StateMap + ?Sized>(
        &self,
        state_map: &S,
        state_buf: &[f64],
        channels: &[OutputChannel],
        out: &mut [f64],
    ) -> Result<(), AnalysisError> {
        let _ = state_map;
        for (slot, &channel) in out.iter_mut().zip(channels) {
            *slot = match channel {
                OutputChannel::State(idx) => *state_buf.get(idx).ok_or(AnalysisError::Invalid {
                    what: "output state channel index out of bounds",
                })?,
                OutputChannel::PortInput { port, component } => {
                    let g = self
                        .input_index(port, component)
                        .ok_or(AnalysisError::Invalid {
                            what: "output PortInput channel out of bounds",
                        })?;
                    self.accepted_inputs[g]
                }
                OutputChannel::PortOutput { port, component } => {
                    let g = self
                        .output_index(port, component)
                        .ok_or(AnalysisError::Invalid {
                            what: "output PortOutput channel out of bounds",
                        })?;
                    self.accepted_outputs[g]
                }
            };
        }
        Ok(())
    }

    /// Linearise the coupled [`System`] about its current accepted state,
    /// returning both the discrete one-step transition `(M_d, B_d, C_d, D_d)`
    /// and the exact Cayley/bilinear continuous-time model `(A, B, C, D)`.
    ///
    /// The system's accepted state, clock, ledger, and boundary biases are
    /// always restored before returning.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any invalid parameter, step, singular Cayley transform, or cancellation error.
    #[allow(clippy::too_many_lines)]
    pub fn linearize_state_space<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        input_channels: &[BoundaryChannel],
        output_channels: &[OutputChannel],
        config: LinearizeConfig,
        cx: &Cx<'_>,
    ) -> Result<StateSpaceLinearization, AnalysisError> {
        self.linearize_state_space_for(
            AnalysisKind::Harmonic,
            state_map,
            input_channels,
            output_channels,
            config,
            cx,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn linearize_state_space_for<S: StateMap + ?Sized>(
        &mut self,
        analysis: AnalysisKind,
        state_map: &S,
        input_channels: &[BoundaryChannel],
        output_channels: &[OutputChannel],
        config: LinearizeConfig,
        cx: &Cx<'_>,
    ) -> Result<StateSpaceLinearization, AnalysisError> {
        self.check_analysis_capabilities(analysis)?;
        let n = state_map.dim();
        let m = input_channels.len();
        let p = output_channels.len();
        if n == 0
            || !config.dt_s.is_finite()
            || config.dt_s <= 0.0
            || !config.fd_eps.is_finite()
            || config.fd_eps <= 0.0
        {
            return Err(AnalysisError::Invalid {
                what: "linearize requires positive state dim, dt_s, and fd_eps",
            });
        }
        let mut resolved_inputs = Vec::with_capacity(m);
        for &ch in input_channels {
            resolved_inputs.push(self.resolve_boundary_channel(ch)?);
        }
        for &ch in output_channels {
            match ch {
                OutputChannel::State(idx) if idx >= n => {
                    return Err(AnalysisError::Invalid {
                        what: "output state channel index out of bounds",
                    });
                }
                OutputChannel::PortInput { port, component }
                    if self.input_index(port, component).is_none() =>
                {
                    return Err(AnalysisError::Invalid {
                        what: "output PortInput channel out of bounds",
                    });
                }
                OutputChannel::PortOutput { port, component }
                    if self.output_index(port, component).is_none() =>
                {
                    return Err(AnalysisError::Invalid {
                        what: "output PortOutput channel out of bounds",
                    });
                }
                _ => {}
            }
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let checkpoint = self.checkpoint()?;
        let saved_bias = self.boundary_bias.clone();
        let mut x0 = vec![0.0_f64; n];
        state_map.read(self, &mut x0)?;

        let result = (|| -> Result<StateSpaceLinearization, AnalysisError> {
            // Base step determines the exact admitted step length `dt_eff`.
            let base_report = self.step(config.dt_s, cx)?;
            let dt_eff = base_report.dt_s;
            let base_boundary: Vec<f64> = resolved_inputs
                .iter()
                .map(|&(b, c)| self.boundary_values[b][c])
                .collect();

            let mut m_d = vec![0.0_f64; n * n];
            let mut b_d = vec![0.0_f64; n * m];
            let mut c_d = vec![0.0_f64; p * n];
            let mut d_d = vec![0.0_f64; p * m];

            let mut x_trial = x0.clone();
            let mut x_plus = vec![0.0_f64; n];
            let mut x_minus = vec![0.0_f64; n];
            let mut y_plus = vec![0.0_f64; p];
            let mut y_minus = vec![0.0_f64; p];

            // Columns of M_d and C_d with respect to state x0[j].
            for j in 0..n {
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                let eps_j = clamp_perturbation(x0[j], config.fd_eps);

                // +eps_j
                self.restore(&checkpoint)?;
                x_trial.copy_from_slice(&x0);
                x_trial[j] = x0[j] + eps_j;
                state_map.write(self, &x_trial)?;
                let rep_p = self.step(dt_eff, cx)?;
                if (rep_p.dt_s - dt_eff).abs() > 1.0e-12 * dt_eff {
                    return Err(AnalysisError::Singular {
                        analysis,
                        what: "step size shortened during state perturbation probe",
                    });
                }
                state_map.read(self, &mut x_plus)?;
                self.read_output_channels(state_map, &x_plus, output_channels, &mut y_plus)?;

                // -eps_j
                self.restore(&checkpoint)?;
                x_trial.copy_from_slice(&x0);
                x_trial[j] = x0[j] - eps_j;
                state_map.write(self, &x_trial)?;
                let rep_m = self.step(dt_eff, cx)?;
                if (rep_m.dt_s - dt_eff).abs() > 1.0e-12 * dt_eff {
                    return Err(AnalysisError::Singular {
                        analysis,
                        what: "step size shortened during state perturbation probe",
                    });
                }
                state_map.read(self, &mut x_minus)?;
                self.read_output_channels(state_map, &x_minus, output_channels, &mut y_minus)?;

                let inv_2eps = 0.5 / eps_j;
                for i in 0..n {
                    m_d[i * n + j] = (x_plus[i] - x_minus[i]) * inv_2eps;
                }
                for i in 0..p {
                    c_d[i * n + j] = (y_plus[i] - y_minus[i]) * inv_2eps;
                }
            }

            // Columns of B_d and D_d with respect to boundary channel k.
            for (k, &(b_idx, comp)) in resolved_inputs.iter().enumerate() {
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                let eps_k = config.fd_eps * base_boundary[k].abs().max(1.0);
                let orig_bias = saved_bias[b_idx][comp];

                // +eps_k
                self.restore(&checkpoint)?;
                self.boundary_bias[b_idx][comp] = orig_bias + eps_k;
                let rep_p = self.step(dt_eff, cx)?;
                if (rep_p.dt_s - dt_eff).abs() > 1.0e-12 * dt_eff {
                    return Err(AnalysisError::Singular {
                        analysis,
                        what: "step size shortened during boundary perturbation probe",
                    });
                }
                state_map.read(self, &mut x_plus)?;
                self.read_output_channels(state_map, &x_plus, output_channels, &mut y_plus)?;

                // -eps_k
                self.restore(&checkpoint)?;
                self.boundary_bias[b_idx][comp] = orig_bias - eps_k;
                let rep_m = self.step(dt_eff, cx)?;
                if (rep_m.dt_s - dt_eff).abs() > 1.0e-12 * dt_eff {
                    return Err(AnalysisError::Singular {
                        analysis,
                        what: "step size shortened during boundary perturbation probe",
                    });
                }
                state_map.read(self, &mut x_minus)?;
                self.read_output_channels(state_map, &x_minus, output_channels, &mut y_minus)?;

                self.boundary_bias[b_idx][comp] = orig_bias;
                let inv_2eps = 0.5 / eps_k;
                for i in 0..n {
                    b_d[i * m + k] = (x_plus[i] - x_minus[i]) * inv_2eps;
                }
                for i in 0..p {
                    d_d[i * m + k] = (y_plus[i] - y_minus[i]) * inv_2eps;
                }
            }

            let (a, b, c, d) = cayley_to_continuous(
                analysis,
                n,
                m,
                p,
                dt_eff,
                output_channels,
                &m_d,
                &b_d,
                &c_d,
                &d_d,
            )?;
            Ok(StateSpaceLinearization {
                state_dim: n,
                input_dim: m,
                output_dim: p,
                dt_s: dt_eff,
                operating_state: x0,
                m_d,
                b_d,
                c_d,
                d_d,
                a,
                b,
                c,
                d,
            })
        })();

        self.boundary_bias = saved_bias;
        let _ = self.restore(&checkpoint);
        result
    }

    /// Compute the linearised continuous `H(i omega)` and discrete
    /// `H_d(e^{i omega dt})` frequency-response transfer matrices about the
    /// current accepted state.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any invalid frequency, singular resolvent, or cancellation error.
    pub fn harmonic_response<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        input_channels: &[BoundaryChannel],
        output_channels: &[OutputChannel],
        frequencies_rad_s: &[f64],
        config: LinearizeConfig,
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::Harmonic)?;
        if frequencies_rad_s.is_empty()
            || !frequencies_rad_s.iter().all(|w| w.is_finite() && *w >= 0.0)
        {
            return Err(AnalysisError::Invalid {
                what: "harmonic_response requires non-empty finite non-negative angular frequencies",
            });
        }
        let lin =
            self.linearize_state_space(state_map, input_channels, output_channels, config, cx)?;
        let n = lin.state_dim;
        let m = lin.input_dim;
        let p = lin.output_dim;

        let mut continuous_transfer = Vec::with_capacity(frequencies_rad_s.len());
        let mut discrete_transfer = Vec::with_capacity(frequencies_rad_s.len());
        let mut mat = vec![C64::ZERO; n * n];
        let mut rhs = vec![C64::ZERO; n];

        for &omega in frequencies_rad_s {
            cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

            // Continuous: H(i omega) = C (i omega I - A)^{-1} B + D.
            for r in 0..n {
                for c_idx in 0..n {
                    let diag = if r == c_idx {
                        C64::new(0.0, omega)
                    } else {
                        C64::ZERO
                    };
                    mat[r * n + c_idx] = diag - C64::from_re(lin.a[r * n + c_idx]);
                }
            }
            let lu_c = lu_complex_with_cancel(&mat, n, || cx.is_cancel_requested()).map_err(
                |err| match err {
                    LuComplexError::Cancelled => AnalysisError::Cancelled,
                    LuComplexError::Singular { .. } => AnalysisError::Singular {
                        analysis: AnalysisKind::Harmonic,
                        what: "continuous resolvent (i omega I - A) is singular at requested frequency",
                    },
                },
            )?;
            let mut h_c = vec![C64::ZERO; p * m];
            for k in 0..m {
                for r in 0..n {
                    rhs[r] = C64::from_re(lin.b[r * m + k]);
                }
                lu_c.solve_with_cancel(&mut rhs, || cx.is_cancel_requested())
                    .map_err(|_| AnalysisError::Cancelled)?;
                for out_i in 0..p {
                    let mut acc = C64::from_re(lin.d[out_i * m + k]);
                    for r in 0..n {
                        acc = acc + C64::from_re(lin.c[out_i * n + r]) * rhs[r];
                    }
                    h_c[out_i * m + k] = acc;
                }
            }
            continuous_transfer.push(h_c);

            // Discrete: H_d(z) = C_d (z I - M_d)^{-1} B_d + D_d, z = exp(i omega dt).
            let theta = omega * lin.dt_s;
            let z = C64::new(theta.cos(), theta.sin());
            for r in 0..n {
                for c_idx in 0..n {
                    let diag = if r == c_idx { z } else { C64::ZERO };
                    mat[r * n + c_idx] = diag - C64::from_re(lin.m_d[r * n + c_idx]);
                }
            }
            let lu_d =
                lu_complex_with_cancel(&mat, n, || cx.is_cancel_requested()).map_err(|err| {
                    match err {
                        LuComplexError::Cancelled => AnalysisError::Cancelled,
                        LuComplexError::Singular { .. } => AnalysisError::Singular {
                            analysis: AnalysisKind::Harmonic,
                            what: "discrete resolvent (exp(i omega dt) I - M_d) is singular",
                        },
                    }
                })?;
            let mut h_d = vec![C64::ZERO; p * m];
            for k in 0..m {
                for r in 0..n {
                    rhs[r] = C64::from_re(lin.b_d[r * m + k]);
                }
                lu_d.solve_with_cancel(&mut rhs, || cx.is_cancel_requested())
                    .map_err(|_| AnalysisError::Cancelled)?;
                for out_i in 0..p {
                    let mut acc = C64::from_re(lin.d_d[out_i * m + k]);
                    for r in 0..n {
                        acc = acc + C64::from_re(lin.c_d[out_i * n + r]) * rhs[r];
                    }
                    h_d[out_i * m + k] = acc;
                }
            }
            discrete_transfer.push(h_d);
        }

        Ok(HarmonicResponse {
            linearization: lin,
            frequencies_rad_s: frequencies_rad_s.to_vec(),
            continuous_transfer,
            discrete_transfer,
        })
    }

    /// Compute the eigenmodes, natural frequencies, and damping ratios of the
    /// linearised continuous-time Jacobian `A` about the current accepted state.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any linearisation, eigensolver, or cancellation refusal.
    pub fn modal_analysis<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        config: ModalConfig,
        cx: &Cx<'_>,
    ) -> Result<ModalAnalysis, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::Modal)?;
        let lin = self.linearize_state_space_for(
            AnalysisKind::Modal,
            state_map,
            &[],
            &[],
            config.linearize,
            cx,
        )?;
        let n = lin.state_dim;
        let a_c64: Vec<C64> = lin.a.iter().map(|&v| C64::from_re(v)).collect();
        let eigs = eig_with_cancel(&a_c64, n, config.max_qr_sweeps, || cx.is_cancel_requested())?;
        let mut modes: Vec<ModalMode> = eigs
            .into_iter()
            .map(|lambda| {
                let omega_d = lambda.im.abs();
                let omega_n = lambda.abs();
                let zeta = if omega_n > 0.0 {
                    -lambda.re / omega_n
                } else {
                    0.0
                };
                ModalMode {
                    eigenvalue: lambda,
                    frequency_hz: omega_d / TAU,
                    angular_frequency_rad_s: omega_d,
                    natural_frequency_rad_s: omega_n,
                    damping_ratio: zeta,
                    decay_rate_per_s: -lambda.re,
                }
            })
            .collect();
        modes.sort_by(|lhs, rhs| {
            lhs.frequency_hz
                .total_cmp(&rhs.frequency_hz)
                .then_with(|| lhs.decay_rate_per_s.total_cmp(&rhs.decay_rate_per_s))
                .then_with(|| lhs.eigenvalue.im.total_cmp(&rhs.eigenvalue.im))
        });
        Ok(ModalAnalysis {
            linearization: lin,
            modes,
        })
    }

    // -----------------------------------------------------------------------
    // 3. Periodic steady state (Poincaré shooting & fs-orbit harmonic balance)
    // -----------------------------------------------------------------------

    /// Find a periodic steady-state orbit `Phi_T(x0*) = x0*` directly over
    /// [`System::step`] via damped Poincaré Newton shooting, returning the
    /// converged trajectory, monodromy matrix `D Phi_T(x0*)`, and Floquet
    /// multipliers.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// convergence, singularity, or cancellation errors.
    #[allow(clippy::too_many_lines)]
    pub fn solve_periodic_shooting<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        initial_guess: Option<&[f64]>,
        config: PeriodicShootingConfig,
        cx: &Cx<'_>,
    ) -> Result<PeriodicShootingSolution, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::PeriodicShooting)?;
        let n = state_map.dim();
        if n == 0
            || !config.period_s.is_finite()
            || config.period_s <= 0.0
            || config.steps_per_period == 0
            || !config.tolerance.is_finite()
            || config.tolerance <= 0.0
            || config.max_iterations == 0
            || !config.fd_eps.is_finite()
            || config.fd_eps <= 0.0
        {
            return Err(AnalysisError::Invalid {
                what: "periodic shooting requires positive dim, period_s, steps_per_period, tolerance, fd_eps, and max_iterations",
            });
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let checkpoint = self.checkpoint()?;
        let mut x0 = vec![0.0_f64; n];
        if let Some(guess) = initial_guess {
            if guess.len() != n || !guess.iter().all(|v| v.is_finite()) {
                return Err(AnalysisError::Invalid {
                    what: "periodic shooting initial_guess must have length dim and finite entries",
                });
            }
            x0.copy_from_slice(guess);
        } else {
            state_map.read(self, &mut x0)?;
        }

        let dt = config.period_s / config.steps_per_period as f64;
        let shoot_once = |sys: &mut System,
                          start: &[f64],
                          record: bool|
         -> Result<(Vec<f64>, Vec<Vec<f64>>), AnalysisError> {
            sys.restore(&checkpoint)?;
            state_map.write(sys, start)?;
            let mut traj = if record {
                let mut v = Vec::with_capacity(config.steps_per_period + 1);
                v.push(start.to_vec());
                v
            } else {
                Vec::new()
            };
            let mut cur = vec![0.0_f64; n];
            let max_substeps_per_bin = 64_usize;
            for _ in 0..config.steps_per_period {
                let mut remaining = dt;
                let mut substeps = 0_usize;
                while remaining > 1.0e-13 * dt {
                    cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                    substeps += 1;
                    if substeps > max_substeps_per_bin {
                        return Err(AnalysisError::NotConverged {
                            analysis: AnalysisKind::PeriodicShooting,
                            iterations: substeps,
                            residual: remaining,
                        });
                    }
                    let rep = sys.step(remaining, cx)?;
                    remaining = (remaining - rep.dt_s).max(0.0);
                }
                if record {
                    state_map.read(sys, &mut cur)?;
                    traj.push(cur.clone());
                }
            }
            state_map.read(sys, &mut cur)?;
            Ok((cur, traj))
        };

        let result = (|| -> Result<PeriodicShootingSolution, AnalysisError> {
            let mut iterations = 0_usize;
            let mut residual = vec![0.0_f64; n];
            let mut monodromy = vec![0.0_f64; n * n];
            let mut jac = vec![0.0_f64; n * n];
            let mut pivots = vec![0_usize; n];
            let mut rhs = vec![0.0_f64; n];
            let mut x_trial = vec![0.0_f64; n];

            let (mut x_end, mut trajectory) = shoot_once(self, &x0, true)?;
            for i in 0..n {
                residual[i] = x_end[i] - x0[i];
            }
            let mut res_norm = l2_norm(&residual);

            while res_norm > config.tolerance {
                if iterations >= config.max_iterations {
                    return Err(AnalysisError::NotConverged {
                        analysis: AnalysisKind::PeriodicShooting,
                        iterations,
                        residual: res_norm,
                    });
                }
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

                // Build monodromy D Phi_T(x0) by central differences.
                for j in 0..n {
                    let eps_j = clamp_perturbation(x0[j], config.fd_eps);
                    x_trial.copy_from_slice(&x0);
                    x_trial[j] = x0[j] + eps_j;
                    let (xp, _) = shoot_once(self, &x_trial, false)?;
                    x_trial[j] = x0[j] - eps_j;
                    let (xm, _) = shoot_once(self, &x_trial, false)?;
                    let inv_2eps = 0.5 / eps_j;
                    for i in 0..n {
                        let m_ij = (xp[i] - xm[i]) * inv_2eps;
                        monodromy[i * n + j] = m_ij;
                        jac[i * n + j] = m_ij - if i == j { 1.0 } else { 0.0 };
                    }
                }
                if !factor_dense(&mut jac, &mut pivots, n) {
                    return Err(AnalysisError::Singular {
                        analysis: AnalysisKind::PeriodicShooting,
                        what: "monodromy Jacobian (D Phi_T - I) is singular",
                    });
                }
                for i in 0..n {
                    rhs[i] = -residual[i];
                }
                if !solve_factored(&jac, &pivots, &mut rhs, n) {
                    return Err(AnalysisError::Singular {
                        analysis: AnalysisKind::PeriodicShooting,
                        what: "Newton update for periodic shooting was non-finite",
                    });
                }

                // Backtracking line search on ||Phi_T(x) - x||_2.
                let mut alpha = 1.0_f64;
                let mut accepted = false;
                for _ in 0..10 {
                    for i in 0..n {
                        x_trial[i] = x0[i] + alpha * rhs[i];
                    }
                    if let Ok((cand_end, cand_traj)) = shoot_once(self, &x_trial, true) {
                        let mut cand_res = vec![0.0_f64; n];
                        for i in 0..n {
                            cand_res[i] = cand_end[i] - x_trial[i];
                        }
                        let cand_norm = l2_norm(&cand_res);
                        if cand_norm.is_finite()
                            && (cand_norm < res_norm || cand_norm <= config.tolerance)
                        {
                            x0.copy_from_slice(&x_trial);
                            x_end = cand_end;
                            trajectory = cand_traj;
                            residual = cand_res;
                            res_norm = cand_norm;
                            accepted = true;
                            break;
                        }
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    return Err(AnalysisError::NotConverged {
                        analysis: AnalysisKind::PeriodicShooting,
                        iterations: iterations + 1,
                        residual: res_norm,
                    });
                }
                iterations += 1;
            }

            // Compute exact monodromy at converged x0* for Floquet multipliers.
            for j in 0..n {
                let eps_j = clamp_perturbation(x0[j], config.fd_eps);
                x_trial.copy_from_slice(&x0);
                x_trial[j] = x0[j] + eps_j;
                let (xp, _) = shoot_once(self, &x_trial, false)?;
                x_trial[j] = x0[j] - eps_j;
                let (xm, _) = shoot_once(self, &x_trial, false)?;
                let inv_2eps = 0.5 / eps_j;
                for i in 0..n {
                    monodromy[i * n + j] = (xp[i] - xm[i]) * inv_2eps;
                }
            }
            let _ = x_end;
            let mono_c64: Vec<C64> = monodromy.iter().map(|&v| C64::from_re(v)).collect();
            let mut floquet = eig_with_cancel(&mono_c64, n, config.max_qr_sweeps, || {
                cx.is_cancel_requested()
            })?;
            floquet.sort_by(|a, b| {
                b.abs()
                    .total_cmp(&a.abs())
                    .then_with(|| a.re.total_cmp(&b.re))
                    .then_with(|| a.im.total_cmp(&b.im))
            });
            let spectral_radius = floquet.iter().map(|mu| mu.abs()).fold(0.0_f64, f64::max);

            Ok(PeriodicShootingSolution {
                initial_state: x0,
                closure_residual: res_norm,
                iterations,
                trajectory,
                monodromy,
                floquet_multipliers: floquet,
                spectral_radius,
            })
        })();

        let _ = self.restore(&checkpoint);
        let sol = result?;
        if config.commit {
            state_map.write(self, &sol.initial_state)?;
        }
        Ok(sol)
    }

    /// Solve for a periodic steady-state orbit directly through `fs-orbit`'s
    /// scaled harmonic-balance engine ([`fs_orbit::solve_hb_scaled_seeded`]),
    /// evaluating the continuous-time vector field `dx/dt = f(t, x)` of this
    /// unchanged [`System`] via `O(dt^2)` Richardson extrapolation over
    /// [`System::step`].
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`,
    /// [`AnalysisError::Orbit`], [`AnalysisError::System`], or
    /// [`AnalysisError::Cancelled`]. Always restores the pre-call checkpoint.
    #[allow(clippy::too_many_lines, clippy::items_after_statements)]
    pub fn solve_periodic_hb<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        config: PeriodicHbConfig,
        cx: &Cx<'_>,
    ) -> Result<PeriodicHbSolution, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::PeriodicHarmonicBalance)?;
        let n = state_map.dim();
        if n == 0 || !config.probe_dt_s.is_finite() || config.probe_dt_s <= 0.0 {
            return Err(AnalysisError::Invalid {
                what: "periodic harmonic balance requires positive state dim and probe_dt_s",
            });
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let checkpoint = self.checkpoint()?;
        let mut x0 = vec![0.0_f64; n];
        state_map.read(self, &mut x0)?;

        let n_harm = config.budget.harmonics;
        let coeff_len = n * (2 * n_harm + 1);
        let needs_omega_slot = !matches!(config.anchor, HbAnchor::Forced { .. });
        let seed = if let Some(s) = config.seed {
            s
        } else {
            let mut s = vec![0.0_f64; coeff_len + usize::from(needs_omega_slot)];
            s[..n].copy_from_slice(&x0);
            if let HbAnchor::Autonomous { omega_guess } | HbAnchor::Backbone { omega_guess, .. } =
                config.anchor
            {
                s[coeff_len] = omega_guess;
            }
            s
        };

        struct SystemOrbitAdapter<'a, 'cx, S: StateMap + ?Sized> {
            system: RefCell<&'a mut System>,
            checkpoint: &'a Checkpoint,
            state_map: &'a S,
            dim: usize,
            probe_dt_s: f64,
            autonomous: bool,
            cx: &'a Cx<'cx>,
            error: RefCell<Option<AnalysisError>>,
        }

        impl<S: StateMap + ?Sized> OrbitProblem for SystemOrbitAdapter<'_, '_, S> {
            fn dim(&self) -> usize {
                self.dim
            }

            fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
                if self.cx.is_cancel_requested() {
                    *self.error.borrow_mut() = Some(AnalysisError::Cancelled);
                    out.fill(f64::NAN);
                    return;
                }
                if self.error.borrow().is_some() {
                    out.fill(f64::NAN);
                    return;
                }
                let eval_rate =
                    |sys: &mut System, dt: f64| -> Result<(Vec<f64>, f64), AnalysisError> {
                        sys.restore(self.checkpoint)?;
                        let mut time_sum = Sum::default();
                        time_sum.add(t);
                        sys.time = time_sum;
                        self.state_map.write(sys, x)?;
                        let rep = sys.step(dt, self.cx)?;
                        let mut x_next = vec![0.0_f64; self.dim];
                        self.state_map.read(sys, &mut x_next)?;
                        for i in 0..self.dim {
                            x_next[i] = (x_next[i] - x[i]) / rep.dt_s;
                        }
                        Ok((x_next, rep.dt_s))
                    };
                let mut sys = self.system.borrow_mut();
                match eval_rate(&mut sys, self.probe_dt_s) {
                    Ok((d_full, dt1)) => match eval_rate(&mut sys, 0.5 * dt1) {
                        Ok((d_half, dt2)) => {
                            let denom = dt1 - dt2;
                            if denom > 0.0 {
                                // General O(dt^2) Richardson extrapolation valid even if `limit_dt_s`
                                // shortened `dt1` or `dt2` unevenly.
                                let w_half = dt1 / denom;
                                let w_full = dt2 / denom;
                                for i in 0..self.dim {
                                    out[i] = w_half * d_half[i] - w_full * d_full[i];
                                }
                            } else {
                                out.copy_from_slice(&d_half);
                            }
                        }
                        Err(e) => {
                            *self.error.borrow_mut() = Some(e);
                            out.fill(f64::NAN);
                        }
                    },
                    Err(e) => {
                        *self.error.borrow_mut() = Some(e);
                        out.fill(f64::NAN);
                    }
                }
            }

            fn autonomous(&self) -> bool {
                self.autonomous
            }
        }

        let adapter = SystemOrbitAdapter {
            system: RefCell::new(self),
            checkpoint: &checkpoint,
            state_map,
            dim: n,
            probe_dt_s: config.probe_dt_s,
            autonomous: matches!(config.anchor, HbAnchor::Autonomous { .. }),
            cx,
            error: RefCell::new(None),
        };
        let hb_res = solve_hb_scaled_seeded(
            &adapter,
            config.anchor,
            seed,
            &config.budget,
            &config.scales,
        );
        let latched = adapter.error.borrow_mut().take();
        drop(adapter);

        let _ = self.restore(&checkpoint);
        if cx.is_cancel_requested() {
            return Err(AnalysisError::Cancelled);
        }
        if let Some(err) = latched {
            return Err(err);
        }
        let orbit = hb_res?;
        Ok(PeriodicHbSolution {
            orbit,
            probe_dt_s: config.probe_dt_s,
        })
    }

    // -----------------------------------------------------------------------
    // 4. Forward tangent & matrix-free IFT adjoint sensitivities
    // -----------------------------------------------------------------------

    /// Compute the exact matrix-free forward tangent `(du, dy)` of one coupled
    /// step along a boundary perturbation `d_boundaries` (in declaration order)
    /// using [`fs_ports::component::PortTangent::jvp`] and `fs-solver`'s GMRES.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any step, dimension, convergence, or cancellation error.
    pub fn step_sensitivity_tangent(
        &mut self,
        dt_s: f64,
        d_boundaries: &[f64],
        config: SensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<StepTangentSensitivity, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::TangentSensitivity)?;
        let nb = self.boundary_scalar_count();
        if d_boundaries.len() != nb
            || !d_boundaries.iter().all(|v| v.is_finite())
            || !config.gmres_tolerance.is_finite()
            || config.gmres_tolerance <= 0.0
            || config.max_gmres_cycles == 0
        {
            return Err(AnalysisError::Invalid {
                what: "step_sensitivity_tangent requires finite d_boundaries matching boundary_scalar_count and positive solver controls",
            });
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let pre_step = self.checkpoint()?;
        let step_res = self.step(dt_s, cx);
        let report = match step_res {
            Ok(rep) => rep,
            Err(err) => {
                let _ = self.restore(&pre_step);
                return Err(err.into());
            }
        };

        let boundary_map = self.flat_boundary_map();
        let n_in = self.inputs.len();
        let mut rhs = vec![0.0_f64; n_in];
        for (k, &(_, _, inp_idx)) in boundary_map.iter().enumerate() {
            // R_i(u; b) = u_i - b_k = 0  ==>  (dR/du) du = -(dR/db) db = +db_k.
            rhs[inp_idx] += d_boundaries[k];
        }

        let op = SystemInterfaceOp::new(
            self,
            &pre_step.accepted,
            report.dt_s,
            &self.accepted_inputs,
            &self.accepted_outputs,
            cx,
        );

        let solve_res = (|| -> Result<StepTangentSensitivity, AnalysisError> {
            let (d_inputs, iters, rel_res) = if n_in == 0 || l2_norm(&rhs) == 0.0 {
                (vec![0.0_f64; n_in], 0, 0.0)
            } else {
                let restart = n_in.clamp(1, 30);
                let mut gmres = GmresState::new(&rhs, restart);
                let rep = gmres.run(
                    &op,
                    &rhs,
                    config.gmres_tolerance,
                    config.max_gmres_cycles,
                    false,
                );
                op.check_status()?;
                if !rep.converged || !gmres.x.iter().all(|v| v.is_finite()) {
                    return Err(AnalysisError::NotConverged {
                        analysis: AnalysisKind::TangentSensitivity,
                        iterations: rep.iters,
                        residual: rep.rel_residual,
                    });
                }
                (gmres.x, rep.iters, rep.rel_residual)
            };
            let mut d_outputs = vec![0.0_f64; self.outputs.len()];
            op.forward_port_jvp(&d_inputs, &mut d_outputs)?;
            Ok(StepTangentSensitivity {
                dt_s: report.dt_s,
                d_inputs,
                d_outputs,
                gmres_iters: iters,
                gmres_residual: rel_res,
            })
        })();

        if !config.commit_step || solve_res.is_err() {
            let _ = self.restore(&pre_step);
        }
        solve_res
    }

    /// Compute the exact matrix-free IFT adjoint gradient `dJ / d b` of one
    /// coupled step directly via [`fs_adjoint::ift_gradient_matfree`] and
    /// [`fs_ports::component::PortTangent::vjp`].
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any invalid input, step, convergence, or cancellation refusal.
    pub fn step_sensitivity_adjoint(
        &mut self,
        dt_s: f64,
        objective: &StepObjective,
        config: SensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<StepAdjointSensitivity, AnalysisError> {
        self.step_sensitivity_adjoint_with_pullback(dt_s, objective, &[], None, config, cx)
    }

    /// General matrix-free IFT adjoint sensitivity of one coupled step via
    /// [`fs_adjoint::ift_gradient_matfree`], additionally pulling back into an
    /// arbitrary caller-defined parameter block (e.g. material-card parameters)
    /// via `parameter_pullback(lambda, bar_y_total) -> dJ_implicit/dp`.
    ///
    /// Specifically, given `R(u; b, p) = u - D(y(u; p), b) = 0`, the closure
    /// receives:
    /// - `lambda`: the solved interface adjoint multiplier `(dR/du)^{-T} (dJ/du + (dy/du)^T dJ/dy)`,
    /// - `bar_y_total = dJ/dy + (dD/dy)^T lambda`: the exact total output cotangent.
    ///
    /// It returns `(dR/dp)^T lambda - (dy/dp)^T dJ/dy` (so `-pullback` equals
    /// `+(dy/dp)^T bar_y_total`, matching `fs_adjoint::ift_gradient_matfree`'s
    /// sign convention `dJ/dp = dj_dp - dr_dp_t(lambda)`!).
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any invalid dimension, step, convergence, or cancellation refusal.
    #[allow(clippy::too_many_lines, clippy::type_complexity)]
    pub fn step_sensitivity_adjoint_with_pullback(
        &mut self,
        dt_s: f64,
        objective: &StepObjective,
        dj_dp: &[f64],
        parameter_pullback: Option<&dyn Fn(&[f64], &[f64]) -> Vec<f64>>,
        config: SensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<StepAdjointSensitivity, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::AdjointSensitivity)?;
        let n_in = self.inputs.len();
        let n_out = self.outputs.len();
        let n_b = self.boundary_scalar_count();
        if (!objective.dj_du.is_empty() && objective.dj_du.len() != n_in)
            || (!objective.dj_dy.is_empty() && objective.dj_dy.len() != n_out)
            || (!objective.dj_db.is_empty() && objective.dj_db.len() != n_b)
            || !objective.dj_du.iter().all(|v| v.is_finite())
            || !objective.dj_dy.iter().all(|v| v.is_finite())
            || !objective.dj_db.iter().all(|v| v.is_finite())
            || !dj_dp.iter().all(|v| v.is_finite())
            || !config.gmres_tolerance.is_finite()
            || config.gmres_tolerance <= 0.0
            || config.max_gmres_cycles == 0
        {
            return Err(AnalysisError::Invalid {
                what: "step_sensitivity_adjoint received mismatched or non-finite objective partials or controls",
            });
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        let pre_step = self.checkpoint()?;
        let report = match self.step(dt_s, cx) {
            Ok(rep) => rep,
            Err(err) => {
                let _ = self.restore(&pre_step);
                return Err(err.into());
            }
        };

        let boundary_map = self.flat_boundary_map();
        let op = SystemInterfaceOp::new(
            self,
            &pre_step.accepted,
            report.dt_s,
            &self.accepted_inputs,
            &self.accepted_outputs,
            cx,
        );

        let solve_res = (|| -> Result<StepAdjointSensitivity, AnalysisError> {
            // Total partial derivative of J(u, y(u)) w.r.t. u:
            // rhs = dj_du + (df/du)^T dj_dy.
            let mut rhs = vec![0.0_f64; n_in];
            if !objective.dj_dy.is_empty() {
                op.reverse_port_vjp(&objective.dj_dy, &mut rhs)?;
            }
            if !objective.dj_du.is_empty() {
                for i in 0..n_in {
                    rhs[i] += objective.dj_du[i];
                }
            }

            // Capture solved lambda inside the pullback passed to `ift_gradient_matfree`
            // so ONE transposed GMRES solve produces both `dJ/db` and `lambda` / `bar_y_total`.
            let captured_lambda: RefCell<Vec<f64>> = RefCell::new(Vec::new());
            let dr_db_t = |lambda: &[f64]| -> Vec<f64> {
                *captured_lambda.borrow_mut() = lambda.to_vec();
                boundary_map
                    .iter()
                    .map(|&(_, _, inp_idx)| -lambda[inp_idx])
                    .collect()
            };

            let (boundary_gradient, adj_report) = if n_in == 0 {
                let g = if objective.dj_db.is_empty() {
                    vec![0.0_f64; n_b]
                } else {
                    objective.dj_db.clone()
                };
                (
                    g,
                    AdjointReport {
                        iters: 0,
                        adjoint_residual: 0.0,
                        converged: true,
                    },
                )
            } else {
                ift_gradient_matfree(
                    &op,
                    &rhs,
                    &objective.dj_db,
                    &dr_db_t,
                    config.gmres_tolerance,
                    config.max_gmres_cycles,
                )
            };

            op.check_status()?;
            if !adj_report.converged || !boundary_gradient.iter().all(|v| v.is_finite()) {
                return Err(AnalysisError::NotConverged {
                    analysis: AnalysisKind::AdjointSensitivity,
                    iterations: adj_report.iters,
                    residual: adj_report.adjoint_residual,
                });
            }

            let adjoint_inputs = captured_lambda.into_inner();
            let mut adjoint_outputs = if objective.dj_dy.is_empty() {
                vec![0.0_f64; n_out]
            } else {
                objective.dj_dy.clone()
            };
            op.scatter_driver_transpose(&adjoint_inputs, &mut adjoint_outputs);

            let parameter_gradient = if let Some(pullback) = parameter_pullback {
                let pb = pullback(&adjoint_inputs, &adjoint_outputs);
                if !dj_dp.is_empty() && dj_dp.len() != pb.len() {
                    return Err(AnalysisError::Invalid {
                        what: "parameter_pullback length did not match dj_dp",
                    });
                }
                let mut g = if dj_dp.is_empty() {
                    vec![0.0_f64; pb.len()]
                } else {
                    dj_dp.to_vec()
                };
                for (gi, pi) in g.iter_mut().zip(&pb) {
                    *gi -= pi;
                }
                if !g.iter().all(|v| v.is_finite()) {
                    return Err(AnalysisError::Singular {
                        analysis: AnalysisKind::AdjointSensitivity,
                        what: "parameter_pullback produced non-finite gradient",
                    });
                }
                g
            } else {
                dj_dp.to_vec()
            };

            Ok(StepAdjointSensitivity {
                dt_s: report.dt_s,
                boundary_gradient,
                parameter_gradient,
                adjoint_inputs,
                adjoint_outputs,
                report: adj_report,
            })
        })();

        if !config.commit_step || solve_res.is_err() {
            let _ = self.restore(&pre_step);
        }
        solve_res
    }

    /// Multi-step reverse-mode adjoint sensitivity over `steps` coupled steps,
    /// using [`System::checkpoint`] / [`System::restore`] across the
    /// two-generation state store and [`fs_adjoint::ift_gradient_matfree`] for
    /// every per-step port pullback.
    ///
    /// Computes exact gradients of a trajectory objective
    /// `J = J_term(x_N) + sum_{k=0}^{N-1} J_k(x_k, u_k, y_k, b_k)` with respect
    /// to the initial accepted state `x_0`, each step's boundary inputs `b_k`,
    /// and a constant boundary shift `b` summed over all steps. Always restores
    /// the initial checkpoint before returning.
    ///
    /// # Errors
    /// [`AnalysisError::MissingTier`] if any component lacks `PortTangent`, or
    /// any invalid input, step, convergence, or cancellation refusal.
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub fn trajectory_adjoint<S: StateMap + ?Sized>(
        &mut self,
        state_map: &S,
        dt_s: f64,
        steps: usize,
        terminal_dj_ds: &[f64],
        step_partials: &[TrajectoryStepPartials],
        config: SensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<TrajectoryAdjointResult, AnalysisError> {
        self.check_analysis_capabilities(AnalysisKind::AdjointSensitivity)?;
        let n = state_map.dim();
        let n_in = self.inputs.len();
        let n_out = self.outputs.len();
        let n_b = self.boundary_scalar_count();
        if n == 0
            || steps == 0
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || terminal_dj_ds.len() != n
            || !terminal_dj_ds.iter().all(|v| v.is_finite())
            || (!step_partials.is_empty() && step_partials.len() != steps)
            || !config.state_fd_eps.is_finite()
            || config.state_fd_eps <= 0.0
        {
            return Err(AnalysisError::Invalid {
                what: "trajectory_adjoint requires positive dim/steps/dt_s, terminal_dj_ds of length dim, and empty or length-steps step_partials",
            });
        }
        for sp in step_partials {
            if (!sp.dj_ds.is_empty() && sp.dj_ds.len() != n)
                || (!sp.dj_du.is_empty() && sp.dj_du.len() != n_in)
                || (!sp.dj_dy.is_empty() && sp.dj_dy.len() != n_out)
                || (!sp.dj_db.is_empty() && sp.dj_db.len() != n_b)
            {
                return Err(AnalysisError::Invalid {
                    what: "trajectory step partials length mismatch",
                });
            }
        }
        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;

        // 1. Forward primal rollout, checkpointing every accepted state.
        let mut checkpoints = Vec::with_capacity(steps);
        let mut step_dts_s = Vec::with_capacity(steps);
        let initial_cp = self.checkpoint()?;
        let saved_bias = self.boundary_bias.clone();

        let rollout_res = (|| -> Result<Vec<f64>, AnalysisError> {
            for _ in 0..steps {
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                checkpoints.push(self.checkpoint()?);
                let rep = self.step(dt_s, cx)?;
                step_dts_s.push(rep.dt_s);
            }
            let mut final_state = vec![0.0_f64; n];
            state_map.read(self, &mut final_state)?;
            Ok(final_state)
        })();
        let final_state = match rollout_res {
            Ok(s) => s,
            Err(err) => {
                self.boundary_bias = saved_bias;
                let _ = self.restore(&initial_cp);
                return Err(err);
            }
        };

        // 2. Reverse sweep over checkpoints[k], k = steps-1 down to 0.
        let boundary_map = self.flat_boundary_map();
        let reverse_res = (|| -> Result<TrajectoryAdjointResult, AnalysisError> {
            let mut bar_s = terminal_dj_ds.to_vec();
            let mut grad_boundary_per_step = vec![vec![0.0_f64; n_b]; steps];
            let mut grad_boundary_total = vec![0.0_f64; n_b];
            let mut step_reports = Vec::with_capacity(steps);

            let mut x_k = vec![0.0_f64; n];
            let mut x_trial = vec![0.0_f64; n];
            let mut x_plus = vec![0.0_f64; n];
            let mut x_minus = vec![0.0_f64; n];

            let mut step_cfg = config;
            step_cfg.commit_step = false;

            for k in (0..steps).rev() {
                cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                let dt_k = step_dts_s[k];
                let sp = step_partials.get(k);

                // (a) Exact matrix-free IFT adjoint for port objective at step k.
                self.restore(&checkpoints[k])?;
                let obj_k = match sp {
                    Some(p) => StepObjective {
                        dj_du: p.dj_du.clone(),
                        dj_dy: p.dj_dy.clone(),
                        dj_db: p.dj_db.clone(),
                    },
                    None => StepObjective {
                        dj_du: Vec::new(),
                        dj_dy: Vec::new(),
                        dj_db: Vec::new(),
                    },
                };
                let step_adj = self.step_sensitivity_adjoint(dt_k, &obj_k, step_cfg, cx)?;
                let mut grad_bk = step_adj.boundary_gradient;
                step_reports.push(step_adj.report);

                // (b) Pull back terminal/downstream state cotangent `bar_s` through
                // the boundary transition `B_k = d x_{k+1} / d b_k` at checkpoint k.
                let bar_s_nonzero = bar_s.iter().any(|&v| v != 0.0);
                if bar_s_nonzero {
                    self.restore(&checkpoints[k])?;
                    self.step(dt_k, cx)?;
                    let base_b: Vec<f64> = boundary_map
                        .iter()
                        .map(|&(b_idx, comp, _)| self.boundary_values[b_idx][comp])
                        .collect();
                    for (b_flat, &(b_idx, comp, _)) in boundary_map.iter().enumerate() {
                        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                        let eps_b = config.state_fd_eps * base_b[b_flat].abs().max(1.0);
                        let orig = saved_bias[b_idx][comp];

                        self.restore(&checkpoints[k])?;
                        self.boundary_bias[b_idx][comp] = orig + eps_b;
                        let rep_p = self.step(dt_k, cx)?;
                        if (rep_p.dt_s - dt_k).abs() > 1.0e-12 * dt_k {
                            return Err(AnalysisError::Singular {
                                analysis: AnalysisKind::AdjointSensitivity,
                                what: "step size shortened during trajectory boundary sensitivity probe",
                            });
                        }
                        state_map.read(self, &mut x_plus)?;

                        self.restore(&checkpoints[k])?;
                        self.boundary_bias[b_idx][comp] = orig - eps_b;
                        let rep_m = self.step(dt_k, cx)?;
                        if (rep_m.dt_s - dt_k).abs() > 1.0e-12 * dt_k {
                            return Err(AnalysisError::Singular {
                                analysis: AnalysisKind::AdjointSensitivity,
                                what: "step size shortened during trajectory boundary sensitivity probe",
                            });
                        }
                        state_map.read(self, &mut x_minus)?;
                        self.boundary_bias[b_idx][comp] = orig;

                        let inv_2eps = 0.5 / eps_b;
                        let mut acc = 0.0_f64;
                        for i in 0..n {
                            acc += bar_s[i] * (x_plus[i] - x_minus[i]) * inv_2eps;
                        }
                        grad_bk[b_flat] += acc;
                    }
                }

                for b_flat in 0..n_b {
                    grad_boundary_total[b_flat] += grad_bk[b_flat];
                }
                grad_boundary_per_step[k] = grad_bk;

                // (c) Pull back `bar_s` and any step-k port objective through `x_k`.
                let has_port_obj = sp.is_some_and(|p| {
                    p.dj_du.iter().any(|&v| v != 0.0) || p.dj_dy.iter().any(|&v| v != 0.0)
                });
                let mut bar_s_prev = vec![0.0_f64; n];
                if let Some(p) = sp
                    && !p.dj_ds.is_empty()
                {
                    bar_s_prev.copy_from_slice(&p.dj_ds);
                }
                if bar_s_nonzero || has_port_obj {
                    self.restore(&checkpoints[k])?;
                    state_map.read(self, &mut x_k)?;
                    let mut u_plus = Vec::new();
                    let mut u_minus = Vec::new();
                    let mut y_plus = Vec::new();
                    let mut y_minus = Vec::new();

                    for j in 0..n {
                        cx.checkpoint().map_err(|_| AnalysisError::Cancelled)?;
                        let eps_j = clamp_perturbation(x_k[j], config.state_fd_eps);

                        self.restore(&checkpoints[k])?;
                        x_trial.copy_from_slice(&x_k);
                        x_trial[j] = x_k[j] + eps_j;
                        state_map.write(self, &x_trial)?;
                        let rep_p = self.step(dt_k, cx)?;
                        if (rep_p.dt_s - dt_k).abs() > 1.0e-12 * dt_k {
                            return Err(AnalysisError::Singular {
                                analysis: AnalysisKind::AdjointSensitivity,
                                what: "step size shortened during trajectory state sensitivity probe",
                            });
                        }
                        state_map.read(self, &mut x_plus)?;
                        if has_port_obj {
                            u_plus.clear();
                            u_plus.extend_from_slice(&self.accepted_inputs);
                            y_plus.clear();
                            y_plus.extend_from_slice(&self.accepted_outputs);
                        }

                        self.restore(&checkpoints[k])?;
                        x_trial.copy_from_slice(&x_k);
                        x_trial[j] = x_k[j] - eps_j;
                        state_map.write(self, &x_trial)?;
                        let rep_m = self.step(dt_k, cx)?;
                        if (rep_m.dt_s - dt_k).abs() > 1.0e-12 * dt_k {
                            return Err(AnalysisError::Singular {
                                analysis: AnalysisKind::AdjointSensitivity,
                                what: "step size shortened during trajectory state sensitivity probe",
                            });
                        }
                        state_map.read(self, &mut x_minus)?;
                        if has_port_obj {
                            u_minus.clear();
                            u_minus.extend_from_slice(&self.accepted_inputs);
                            y_minus.clear();
                            y_minus.extend_from_slice(&self.accepted_outputs);
                        }

                        let inv_2eps = 0.5 / eps_j;
                        let mut acc = 0.0_f64;
                        if bar_s_nonzero {
                            for i in 0..n {
                                acc += bar_s[i] * (x_plus[i] - x_minus[i]) * inv_2eps;
                            }
                        }
                        if let Some(p) = sp {
                            if !p.dj_du.is_empty() {
                                for i in 0..n_in {
                                    acc += p.dj_du[i] * (u_plus[i] - u_minus[i]) * inv_2eps;
                                }
                            }
                            if !p.dj_dy.is_empty() {
                                for i in 0..n_out {
                                    acc += p.dj_dy[i] * (y_plus[i] - y_minus[i]) * inv_2eps;
                                }
                            }
                        }
                        bar_s_prev[j] += acc;
                    }
                }
                bar_s = bar_s_prev;
            }
            step_reports.reverse();

            Ok(TrajectoryAdjointResult {
                steps,
                step_dts_s,
                final_state,
                grad_initial_state: bar_s,
                grad_boundary_per_step,
                grad_boundary_total,
                step_reports,
            })
        })();

        self.boundary_bias = saved_bias;
        let _ = self.restore(&initial_cp);
        reverse_res
    }
}
