//! Discrete adjoints of the accepted reference-mass enthalpy balance.
//!
//! For `R=M(h-h_old)+dt*(A(T(h),f(h))*T(h)-b(q))`, solve
//! `J_h^T lambda=h_bar`, then return `M*lambda` for history and
//! `dt*M_source^T*lambda` for nodal volumetric source density. Here
//! the spatial `J_h` already includes temperature-column scaling and any
//! declared phase-conductivity chain rule. The shared bounded FGMRES solver
//! sees its explicit transpose action and the exact diagonal Jacobi inverse.
//! The ambient-radiation binder can retain physical rank-one feedback terms;
//! these participate in both actions, Jacobi and all history/source pullbacks.
//!
//! Geometry, reference masses, chart, conductivity laws, interfaces and time
//! step are frozen. The base pullback also fixes boundary data; `robin_response`
//! explicitly selects convective references and coefficients. Source derivatives describe the P1
//! nodal source field, including a uniform source represented by equal nodal
//! values. No derivative of solver iterations, phase-law parameters, material
//! parameters, moving geometry or chart selection is inferred. Chart corners
//! and validity endpoints refuse classical two-sided derivatives; a latent
//! plateau interior has exactly zero temperature sensitivity.

use std::{cell::RefCell, fmt};

use fs_exec::Cx;
use fs_solver::{
    FgmresState, FlexiblePreconditioner, LinearOp, NewtonKrylovState, SolverRunProgress, norm2,
};
use fs_sparse::Csr;

/// Convection-only response of the checked total-enthalpy endpoint.
pub mod robin;
pub use robin::{EnthalpyRobinGradient, EnthalpyRobinResponse};

use super::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
    poll, validate_work,
};
use crate::{
    ConductionError, ConductionMesh, ConductionProblem, LinearConfig, ThermalBc, ThermalInterfaces,
    assemble::ASSEMBLY_TILE,
};

/// A refused derivative publishes no partial gradient or accepted replacement.
#[derive(Debug, Clone, PartialEq)]
pub enum EnthalpyAdjointError {
    /// Shape, solve budget, or diagonal admission failed.
    InvalidInput(&'static str),
    /// Primal chart, spatial assembly, energy, or cancellation refusal.
    Enthalpy(EnthalpyError),
    /// The supplied endpoint fails the recomputed original Newton target.
    PrimalResidual {
        /// Actual Euclidean residual norm, J.
        residual_j: f64,
        /// Original absolute/relative Newton target, J.
        tolerance_j: f64,
    },
    /// Adjacent chart segments have different temperature slopes or different
    /// fraction slopes used by an active phase-conductivity law.
    ChartKink {
        /// Vertex on the ambiguous branch boundary.
        vertex: usize,
        /// Specific enthalpy at that boundary, J/kg.
        specific_enthalpy_j_kg: f64,
    },
    /// A finite chart endpoint has no two-sided constitutive neighborhood.
    ChartEndpoint {
        /// Vertex at the edge of the declared chart.
        vertex: usize,
        /// Specific enthalpy at the chart endpoint, J/kg.
        specific_enthalpy_j_kg: f64,
    },
    /// A sampled conductivity law has a slope corner or validity endpoint.
    MaterialKink,
    /// The bounded transpose solve failed its independently recomputed gate.
    NotConverged {
        /// Actual inner Krylov columns.
        iterations: usize,
        /// True Euclidean relative residual.
        relative_residual: f64,
        /// Caller-requested relative tolerance.
        tolerance: f64,
    },
    /// A finite input produced an unrepresentable derivative quantity.
    NonFiniteArithmetic,
}

impl fmt::Display for EnthalpyAdjointError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(reason) => write!(f, "enthalpy adjoint input: {reason}"),
            Self::Enthalpy(error) => error.fmt(f),
            Self::PrimalResidual {
                residual_j,
                tolerance_j,
            } => write!(
                f,
                "enthalpy endpoint residual {residual_j} J exceeds {tolerance_j} J"
            ),
            Self::ChartKink {
                vertex,
                specific_enthalpy_j_kg,
            } => write!(
                f,
                "enthalpy derivative at vertex {vertex}, h={specific_enthalpy_j_kg}, crosses a chart slope corner"
            ),
            Self::ChartEndpoint {
                vertex,
                specific_enthalpy_j_kg,
            } => write!(
                f,
                "enthalpy derivative at vertex {vertex}, h={specific_enthalpy_j_kg}, lies at a chart validity endpoint"
            ),
            Self::MaterialKink => {
                write!(f, "conductivity derivative has no unique two-sided slope")
            }
            Self::NotConverged {
                iterations,
                relative_residual,
                tolerance,
            } => write!(
                f,
                "enthalpy adjoint stopped after {iterations} columns: residual {relative_residual}, tolerance {tolerance}"
            ),
            Self::NonFiniteArithmetic => {
                write!(f, "enthalpy derivative arithmetic is not representable")
            }
        }
    }
}

impl std::error::Error for EnthalpyAdjointError {}

impl From<EnthalpyError> for EnthalpyAdjointError {
    fn from(error: EnthalpyError) -> Self {
        Self::Enthalpy(error)
    }
}

impl From<ConductionError> for EnthalpyAdjointError {
    fn from(error: ConductionError) -> Self {
        Self::Enthalpy(error.into())
    }
}

/// Pullbacks of a scalar endpoint objective, plus the actual transpose solve.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpyStepGradient {
    /// Derivative with respect to each previous nodal specific enthalpy.
    pub previous_specific_enthalpy: Vec<f64>,
    /// Derivative with respect to each P1 nodal source density [W/m3].
    /// Sum these entries for a spatially uniform source-amplitude derivative.
    pub source_density: Vec<f64>,
    /// Multiplier of the physical residual in joules, `J_h^T lambda=h_bar`.
    pub adjoint: Vec<f64>,
    /// Independently recomputed true Euclidean relative transpose residual.
    pub relative_residual: f64,
    /// Actual inner Krylov columns, excluding the primal solve.
    pub iterations: usize,
}

/// An accepted, rechecked endpoint with an owned sparse physical tangent and
/// optional retained mean-radiation feedback factors.
/// Only the immutable mesh is borrowed; step-local source/boundary fields may
/// be dropped. A caller can retain a bounded vector of these for a small
/// full-storage reverse sweep without retaining any Newton/Krylov iteration.
#[derive(Debug)]
pub struct EnthalpyStepLinearization<'m> {
    primal: EnthalpyStepSolution,
    mesh: &'m ConductionMesh,
    tangent: Csr,
    slopes: Vec<f64>,
    masses: Vec<f64>,
    inverse_diagonal: Vec<f64>,
    dt: f64,
    // Only the physical ambient-radiation binder installs these updates.
    // They are in residual-joule/enthalpy coordinates, after chart scaling.
    feedback: Vec<EnthalpyFeedback>,
    // Retained by the checked producer. Radiation restores its ORIGINAL
    // convection boundary after preparing the complete physical tangent.
    convection_boundary: crate::ThermalBoundary,
}

#[derive(Debug)]
struct EnthalpyFeedback {
    left: Vec<f64>,
    right: Vec<f64>,
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Solve the primal and bind its implicit derivatives. The accepted state
    /// is rechecked with the same residual/energy gates as `linearize_accepted`.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let accepted = self.advance(cx, problem, interfaces, old_h, dt_s, config)?;
        self.linearize_accepted(cx, problem, interfaces, old_h, dt_s, config, accepted)
    }

    /// Recompute the physical residual, energy balance and chart fields before
    /// attaching a tangent to a supplied endpoint. Public solution fields and
    /// convergence flags are not trusted as acceptance evidence. The original
    /// Newton target is reconstructed from the residual at `old_h`; all shared
    /// Newton controls are admitted by the shared constructor, without running
    /// another Newton iteration. Returned telemetry retains the supplied
    /// iteration history, with its endpoint residual/convergence recomputed.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn linearize_accepted(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        accepted: EnthalpyStepSolution,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        self.linearize_accepted_with_target(
            cx, problem, interfaces, old_h, dt_s, config, accepted, None,
        )
    }

    // The radiative binder supplies the original complete nonlinear initial
    // residual target, not a target derived from a frozen endpoint secant.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn linearize_accepted_with_target(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        accepted: EnthalpyStepSolution,
        physical_target: Option<f64>,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        let n = self.masses.len();
        if !std::ptr::eq(self.mesh, problem.mesh)
            || old_h.len() != n
            || accepted.specific_enthalpy_j_kg.len() != n
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "history, endpoint and reference mass must match the exact mesh",
            ));
        }
        if problem
            .boundary
            .conditions()
            .iter()
            .any(|bc| matches!(bc, ThermalBc::Dirichlet { .. }))
        {
            return Err(EnthalpyError::UnsupportedDirichlet.into());
        }
        validate_work(n, dt_s, config)?;
        for condition in problem.boundary.conditions() {
            condition.validate(n)?;
        }
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
        };
        let stage = context.stage(&accepted.specific_enthalpy_j_kg)?;
        let initial = NewtonKrylovState::new(&stage, old_h.to_vec(), config.newton)
            .map_err(|error| stage.take_failure().unwrap_or(EnthalpyError::Newton(error)))?;
        let tolerance_j = finite(physical_target.unwrap_or_else(|| {
            config.newton.absolute_tolerance.max(
                config.newton.relative_tolerance * initial.residual_norm().max(f64::MIN_POSITIVE),
            )
        }))?;
        if tolerance_j < 0.0 {
            return Err(EnthalpyAdjointError::InvalidInput(
                "negative physical residual target",
            ));
        }
        drop(initial);
        let mut residual = vec![0.0; n];
        context.residual(&accepted.specific_enthalpy_j_kg, &mut residual)?;
        let residual_j = finite(norm2(&residual))?;
        if residual_j > tolerance_j {
            return Err(EnthalpyAdjointError::PrimalResidual {
                residual_j,
                tolerance_j,
            });
        }
        let mut report = accepted.newton;
        report.residual_norm = residual_j;
        report.converged = true;
        report.diagnosis = None;
        let primal = context.finish(
            accepted.specific_enthalpy_j_kg,
            report,
            config.energy_tolerance_j,
        )?;
        check_chart(cx, self, &primal.specific_enthalpy_j_kg)?;
        if !crate::adjoint::robin::nonlinear::material_is_smooth(cx, problem, &primal.temperature)?
        {
            return Err(EnthalpyAdjointError::MaterialKink);
        }
        let mut inverse_diagonal = Vec::with_capacity(n);
        for row in 0..n {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            let diagonal = finite(self.masses[row] + dt_s * stage.tangent.get(row, row))?;
            if diagonal == 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "nonzero exact enthalpy Jacobian diagonal required for Jacobi",
                ));
            }
            inverse_diagonal.push(finite(1.0 / diagonal)?);
        }
        poll(cx, n)?;
        Ok(EnthalpyStepLinearization {
            primal,
            mesh: self.mesh,
            tangent: stage.tangent,
            slopes: stage.slope,
            masses: self.masses.clone(),
            inverse_diagonal,
            dt: dt_s,
            feedback: Vec::new(),
            convection_boundary: problem.boundary.clone(),
        })
    }
}

impl EnthalpyStepLinearization<'_> {
    pub(super) fn retain_convection_boundary(&mut self, boundary: &crate::ThermalBoundary) {
        self.convection_boundary = boundary.clone();
    }

    /// Attach one physical boundary feedback term. Entry admission is owned by
    /// the radiative binder; no external caller may modify the checked tangent.
    pub(super) fn add_radiation_feedback(
        &mut self,
        cx: &Cx<'_>,
        left: Vec<f64>,
        right: Vec<f64>,
    ) -> Result<(), EnthalpyAdjointError> {
        let n = self.masses.len();
        vector(cx, &left, n)?;
        vector(cx, &right, n)?;
        for i in 0..n {
            if i % ASSEMBLY_TILE == 0 {
                poll(cx, i)?;
            }
            let diagonal = finite(1.0 / self.inverse_diagonal[i] + left[i] * right[i])?;
            if diagonal == 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "nonzero radiative enthalpy Jacobian diagonal required for Jacobi",
                ));
            }
            self.inverse_diagonal[i] = finite(1.0 / diagonal)?;
        }
        self.feedback.push(EnthalpyFeedback { left, right });
        poll(cx, n)?;
        Ok(())
    }

    /// Physically rechecked endpoint and independently integrated energy audit.
    #[must_use]
    pub const fn primal(&self) -> &EnthalpyStepSolution {
        &self.primal
    }

    /// Convert a nodal temperature cotangent to a specific-enthalpy cotangent.
    /// Add any direct enthalpy objective/carry cotangent before `pullback`.
    pub fn temperature_pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        vector(cx, seed, self.masses.len())?;
        let mut result = Vec::with_capacity(seed.len());
        for (vertex, (&bar, &slope)) in seed.iter().zip(&self.slopes).enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            result.push(finite(bar * slope)?);
        }
        poll(cx, seed.len())?;
        Ok(result)
    }

    /// Apply the exact endpoint Jacobian in residual-joule coordinates.
    pub fn apply_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.action(cx, direction, false)
    }

    /// Apply its sparse/feedback transpose, including the correct slope side.
    pub fn apply_jacobian_transpose(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.action(cx, direction, true)
    }

    fn action(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        transpose: bool,
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        vector(cx, direction, self.masses.len())?;
        let mut output = vec![0.0; direction.len()];
        self.action_into(cx, direction, &mut output, transpose)?;
        Ok(output)
    }

    fn action_into(
        &self,
        cx: &Cx<'_>,
        x: &[f64],
        y: &mut [f64],
        transpose: bool,
    ) -> Result<(), EnthalpyAdjointError> {
        y.fill(0.0);
        for (row, &x_row) in x.iter().enumerate() {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            let (columns, entries) = self.tangent.row(row);
            if transpose {
                for (&column, &entry) in columns.iter().zip(entries) {
                    y[column] = finite((entry * x_row).mul_add(self.dt, y[column]))?;
                }
            } else {
                let mut sum = 0.0;
                for (&column, &entry) in columns.iter().zip(entries) {
                    sum = finite(entry.mul_add(x[column], sum))?;
                }
                y[row] = finite(self.dt * sum)?;
            }
        }
        for (row, value) in y.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            *value = finite(self.masses[row].mul_add(x[row], *value))?;
        }
        for update in &self.feedback {
            let (left, right) = if transpose {
                (&update.right, &update.left)
            } else {
                (&update.left, &update.right)
            };
            let mut contraction = 0.0;
            for (i, (&weight, &value)) in right.iter().zip(x).enumerate() {
                if i % ASSEMBLY_TILE == 0 {
                    poll(cx, i)?;
                }
                contraction = finite(weight.mul_add(value, contraction))?;
            }
            for (i, (value, &weight)) in y.iter_mut().zip(left).enumerate() {
                if i % ASSEMBLY_TILE == 0 {
                    poll(cx, i)?;
                }
                *value = finite(weight.mul_add(contraction, *value))?;
            }
        }
        poll(cx, x.len())?;
        Ok(())
    }

    /// Solve the implicit discrete transpose and return history/source VJPs.
    /// `max_iterations` is an exact inner-column cap here, not a cycle count.
    /// Restart is capped by dimension and remaining work. Checked basis and
    /// Hessenberg arithmetic precedes any Krylov allocation. Cancellation is
    /// polled inside both sparse actions and Jacobi, and between cycles.
    #[allow(clippy::too_many_lines)]
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        config: LinearConfig,
    ) -> Result<EnthalpyStepGradient, EnthalpyAdjointError> {
        let n = self.masses.len();
        vector(cx, seed, n)?;
        admit_linear(n, config)?;
        let scale = seed.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let normalized: Vec<f64> = if scale == 0.0 {
            vec![0.0; n]
        } else {
            seed.iter().map(|x| x / scale).collect()
        };
        let failure = RefCell::new(None);
        let operator = Transpose {
            linearization: self,
            cx,
            failure: &failure,
        };
        let restart = config.restart.min(n).min(config.max_iterations);
        let mut state = FgmresState::new(&normalized, restart);
        if scale != 0.0 {
            while state.rel_residual() >= config.tolerance && state.iters < config.max_iterations {
                poll(cx, state.iters)?;
                let before = state.iters;
                state.restart = restart.min(config.max_iterations - before);
                let progress = state.run_cancellable(
                    &operator,
                    &operator,
                    &normalized,
                    config.tolerance,
                    1,
                    cx,
                );
                match progress.progress {
                    SolverRunProgress::Complete => {}
                    SolverRunProgress::Paused => {
                        return Err(ConductionError::Cancelled {
                            stage: "spatial-enthalpy-adjoint",
                            at: state.iters,
                        }
                        .into());
                    }
                    SolverRunProgress::CallbackPanicked(callback) => {
                        return Err(EnthalpyError::SolverCallbackPanicked(callback).into());
                    }
                    SolverRunProgress::DimensionMismatch { .. } => {
                        return Err(EnthalpyAdjointError::InvalidInput(
                            "enthalpy transpose dimension changed",
                        ));
                    }
                }
                if let Some(error) = failure.take() {
                    return Err(error);
                }
                if state.iters == before {
                    break;
                }
            }
        }
        poll(cx, state.iters)?;
        let mut residual = self.apply_jacobian_transpose(cx, &state.x)?;
        for (value, &rhs) in residual.iter_mut().zip(&normalized) {
            *value = finite(rhs - *value)?;
        }
        let relative_residual =
            finite(norm2(&residual) / norm2(&normalized).max(f64::MIN_POSITIVE))?;
        if relative_residual >= config.tolerance {
            return Err(EnthalpyAdjointError::NotConverged {
                iterations: state.iters,
                relative_residual,
                tolerance: config.tolerance,
            });
        }
        let mut previous_specific_enthalpy = Vec::with_capacity(n);
        for (row, value) in state.x.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 {
                poll(cx, row)?;
            }
            *value = finite(*value * scale)?;
            previous_specific_enthalpy.push(finite(self.masses[row] * *value)?);
        }
        let mut source_density = vec![0.0; n];
        for element in 0..self.mesh.element_count() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let tet = self.mesh.complex().tets[element].map(|vertex| vertex as usize);
            let volume = self.mesh.element_volume(element);
            // Same consistent P1 source mass as the existing source assembler;
            // storage masses are lumped, but the source derivative is not.
            for (b, &vertex) in tet.iter().enumerate() {
                let mut value = 0.0;
                for (a, &row) in tet.iter().enumerate() {
                    let mass = if a == b { volume / 10.0 } else { volume / 20.0 };
                    value = finite(mass.mul_add(state.x[row], value))?;
                }
                source_density[vertex] = finite(self.dt.mul_add(value, source_density[vertex]))?;
            }
        }
        poll(cx, n)?;
        Ok(EnthalpyStepGradient {
            previous_specific_enthalpy,
            source_density,
            adjoint: state.x,
            relative_residual,
            iterations: state.iters,
        })
    }
}

struct Transpose<'a, 'm, 'cx> {
    linearization: &'a EnthalpyStepLinearization<'m>,
    cx: &'a Cx<'cx>,
    failure: &'a RefCell<Option<EnthalpyAdjointError>>,
}

impl Transpose<'_, '_, '_> {
    fn retain(&self, result: Result<(), EnthalpyAdjointError>, output: &mut [f64]) {
        if let Err(error) = result {
            if self.failure.borrow().is_none() {
                *self.failure.borrow_mut() = Some(error);
            }
            output.fill(f64::NAN);
        }
    }
}

impl LinearOp for Transpose<'_, '_, '_> {
    fn n(&self) -> usize {
        self.linearization.masses.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.retain(self.linearization.action_into(self.cx, x, y, true), y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.retain(self.linearization.action_into(self.cx, x, y, false), y);
    }
}

impl FlexiblePreconditioner for Transpose<'_, '_, '_> {
    fn apply(&self, _iteration: usize, residual: &[f64], output: &mut [f64]) {
        let result = (|| {
            for (row, value) in output.iter_mut().enumerate() {
                if row % ASSEMBLY_TILE == 0 {
                    poll(self.cx, row)?;
                }
                *value = finite(residual[row] * self.linearization.inverse_diagonal[row])?;
            }
            Ok(())
        })();
        self.retain(result, output);
    }
}

// Exact segment equality is intentional: a tolerance would silently admit a
// genuine corner and claim a classical derivative for a nonsmooth chart.
#[allow(clippy::float_cmp)]
fn check_chart(
    cx: &Cx<'_>,
    storage: &EnthalpyBackwardEuler<'_, '_>,
    h: &[f64],
) -> Result<(), EnthalpyAdjointError> {
    for (vertex, &value) in h.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(cx, vertex)?;
        }
        let curve = storage.curve_for_vertex(vertex);
        let knots = curve.knots();
        if let Ok(index) = knots.binary_search_by(|knot| {
            knot.specific_enthalpy_j_kg
                .partial_cmp(&value)
                .expect("finite admitted chart state")
        }) {
            if index == 0 || index + 1 == knots.len() {
                return Err(EnthalpyAdjointError::ChartEndpoint {
                    vertex,
                    specific_enthalpy_j_kg: value,
                });
            }
            let slope = |h| {
                curve
                    .temperature_derivative_at_specific_enthalpy(h)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })
            };
            let temperature_kink = slope(knots[index - 1].specific_enthalpy_j_kg)? != slope(value)?;
            let mut fraction_kink = false;
            if storage
                .phase_conductivity
                .as_ref()
                .is_some_and(|policy| policy.variable_vertices[vertex])
            {
                let fraction_slope = |h| {
                    curve
                        .liquid_mass_fraction_derivative_at_specific_enthalpy(h)
                        .map_err(|source| EnthalpyError::Phase { vertex, source })
                };
                fraction_kink = fraction_slope(knots[index - 1].specific_enthalpy_j_kg)?
                    != fraction_slope(value)?;
            }
            if temperature_kink || fraction_kink {
                return Err(EnthalpyAdjointError::ChartKink {
                    vertex,
                    specific_enthalpy_j_kg: value,
                });
            }
        }
    }
    Ok(())
}

fn vector(cx: &Cx<'_>, values: &[f64], n: usize) -> Result<(), EnthalpyAdjointError> {
    poll(cx, 0)?;
    if values.len() != n {
        return Err(EnthalpyAdjointError::InvalidInput(
            "nodal vector length mismatch",
        ));
    }
    for (index, &value) in values.iter().enumerate() {
        if index % ASSEMBLY_TILE == 0 {
            poll(cx, index)?;
        }
        finite(value)?;
    }
    Ok(())
}

fn admit_linear(n: usize, config: LinearConfig) -> Result<(), EnthalpyAdjointError> {
    if config.restart == 0
        || config.max_iterations == 0
        || !config.tolerance.is_finite()
        || config.tolerance <= 0.0
        || config.tolerance >= 1.0
    {
        return Err(EnthalpyAdjointError::InvalidInput(
            "positive restart/inner-column budgets and tolerance in (0,1) required",
        ));
    }
    let r = config.restart;
    r.checked_add(1)
        .and_then(|v| v.checked_mul(r))
        .and_then(|h| {
            n.checked_mul(r.checked_mul(2)?.checked_add(16)?)?
                .checked_add(h)
        })
        .and_then(|v| v.checked_add(r.checked_mul(8)?))
        .and_then(|v| v.checked_add(config.max_iterations))
        .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or(EnthalpyAdjointError::InvalidInput(
            "adjoint workspace size overflow",
        ))?;
    Ok(())
}

fn finite(value: f64) -> Result<f64, EnthalpyAdjointError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyAdjointError::NonFiniteArithmetic)
    }
}
