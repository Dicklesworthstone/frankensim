//! Sparse nonlinear heat transport on the original fixed thermal graph.
use super::*;
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolverRunProgress};
use std::cell::RefCell;

/// Bounded linear work inside the original Newton/time-step controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionLinearControls {
    /// Arnoldi dimension, 1..=128. Storage is O(cells * restart + links).
    pub restart: usize,
    /// Restart cycles allowed for one Newton direction before interval retry.
    pub maximum_cycles_per_solve: usize,
    /// Whole-advance ceiling, including rejected nonlinear and adaptive trials.
    pub maximum_iterations: usize,
    /// True Euclidean relative linear residual, strictly between zero and one.
    /// This never replaces the original local and total energy acceptance gates.
    pub relative_tolerance: f64,
}
impl Default for ImplicitConductionLinearControls {
    fn default() -> Self {
        Self {
            restart: 32,
            maximum_cycles_per_solve: 32,
            maximum_iterations: 1_048_576,
            relative_tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) enum LinearTopology {
    Chain,
    Graph(ImplicitConductionLinearControls),
}

/// Branched, cyclic or disconnected fixed-mass thermal networks. All links and
/// reservoirs must have original geometric/material conduction bindings.
/// Parallel interfaces and asymmetric covered areas retain their original
/// meaning. No connectivity or cell ordering is inferred from numeric IDs.
///
/// This uses the same enthalpy residual, analytic conductivity/temperature
/// derivatives, line search, events, heat/entropy accounts and adaptive owner
/// as `ImplicitConductionChain`. Only its linear algebra differs: bounded,
/// diagonally preconditioned FGMRES replaces the chain's O(n) direct solve.
pub struct ImplicitConductionGraph<'a> {
    inner: ImplicitConductionChain<'a>,
}
impl<'a> ImplicitConductionGraph<'a> {
    /// Bind fixed cells/faces and the original caloric/conductivity objects.
    /// Invalid work controls, source domains or incomplete bindings refuse
    /// before any physical state changes.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        conduction: &ConductionObserver<'a>,
        laws: &[(MaterialId, &'a dyn ImplicitConductionLaw)],
        linear: ImplicitConductionLinearControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        if linear.restart == 0
            || linear.restart > 128
            || linear.maximum_cycles_per_solve == 0
            || linear.maximum_iterations == 0
            || !positive(linear.relative_tolerance)
            || linear.relative_tolerance >= 1.0
        {
            return Err(ImplicitConductionError::Invalid(
                "invalid conduction Krylov controls",
            ));
        }
        Ok(Self {
            inner: ImplicitConductionChain::prepare(
                thermal,
                conduction,
                laws,
                LinearTopology::Graph(linear),
                checkpoint,
            )?,
        })
    }

    /// Stable energy-unknown ordering, independent of graph connectivity.
    pub fn cell_ids(&self) -> &[CellId] {
        self.inner.cell_ids()
    }

    /// Advance to absolute seconds using the original nonlinear endpoint solve.
    /// Time accuracy remains caller-controlled; failure returns no candidate.
    pub fn advance_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
        self.inner
            .advance_to(thermal, target_s, controls, workspace, checkpoint)
    }

    /// Advance with the original enthalpy/temperature/phase step-doubling owner.
    /// All coarse, fine and rejected trials share one Krylov iteration budget.
    pub fn advance_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: adaptive::AdaptiveImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            adaptive::AdaptiveImplicitConductionReport,
        ),
        adaptive::AdaptiveImplicitConductionError,
    > {
        self.inner
            .advance_adaptive_to(thermal, target_s, controls, workspace, checkpoint)
    }
}

// The physical Jacobian has one diagonal and two off-diagonals per link, even
// with unequal masses, nonlinear k or latent plateaus (generally nonsymmetric).
// No dense matrix, symmetry assumption or numerical constitutive difference.
struct GraphOperator<'a, 'c> {
    diagonal: &'a [f64],
    entries: &'a [(usize, usize, f64)],
    checkpoint: &'a RefCell<&'c mut dyn FnMut() -> Result<(), DepositionError>>,
    error: &'a RefCell<Option<DepositionError>>,
}
impl GraphOperator<'_, '_> {
    fn poll(&self) -> bool {
        if self.error.borrow().is_some() {
            return false;
        }
        if let Err(error) = (self.checkpoint.borrow_mut())() {
            *self.error.borrow_mut() = Some(error);
            return false;
        }
        true
    }
    fn finite(&self, value: f64) -> bool {
        if value.is_finite() {
            return true;
        }
        *self.error.borrow_mut() = Some(DepositionError::Arithmetic);
        false
    }
    fn action(&self, x: &[f64], y: &mut [f64], transposed: bool) {
        for (i, ((out, &value), &diagonal)) in y.iter_mut().zip(x).zip(self.diagonal).enumerate() {
            if i % 256 == 0 && !self.poll() {
                return;
            }
            *out = diagonal * value;
            if !self.finite(*out) {
                return;
            }
        }
        for (i, &(row, col, value)) in self.entries.iter().enumerate() {
            if i % 256 == 0 && !self.poll() {
                return;
            }
            let (row, col) = if transposed { (col, row) } else { (row, col) };
            y[row] += value * x[col];
            if !self.finite(y[row]) {
                return;
            }
        }
    }
}
impl LinearOp for GraphOperator<'_, '_> {
    fn n(&self) -> usize {
        self.diagonal.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, false);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, true);
    }
}
impl FlexiblePreconditioner for GraphOperator<'_, '_> {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        for (i, ((out, &r), &d)) in output
            .iter_mut()
            .zip(residual)
            .zip(self.diagonal)
            .enumerate()
        {
            if i % 256 == 0 && !self.poll() {
                return;
            }
            // A zero diagonal does not prove a nonsymmetric matrix singular.
            // Leave that coordinate unscaled; FGMRES and the true residual
            // decide whether this original system can be solved.
            *out = if d.abs() >= f64::MIN_POSITIVE {
                r / d
            } else {
                r
            };
            if !self.finite(*out) {
                return;
            }
        }
    }
}

pub(super) fn solve(
    workspace: &mut ImplicitConductionWorkspace,
    controls: ImplicitConductionLinearControls,
    report: &mut ImplicitConductionReport,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ImplicitConductionError> {
    let error = RefCell::new(None);
    let checkpoint: &mut dyn FnMut() -> Result<(), DepositionError> = checkpoint;
    let checkpoint = RefCell::new(checkpoint);
    let operator = GraphOperator {
        diagonal: &workspace.diagonal,
        entries: workspace
            .graph_entries
            .as_deref()
            .ok_or(ImplicitConductionError::Invalid("missing graph Jacobian"))?,
        checkpoint: &checkpoint,
        error: &error,
    };
    let mut linear = FgmresState::new(&workspace.rhs, controls.restart.min(operator.n()));
    for _ in 0..controls.maximum_cycles_per_solve {
        let remaining = controls
            .maximum_iterations
            .saturating_sub(report.linear_iterations);
        if remaining == 0 {
            return Err(DepositionError::Budget("conduction Krylov iterations").into());
        }
        linear.restart = controls.restart.min(operator.n()).min(remaining);
        let before = linear.iters;
        let outcome = linear.run_controlled(
            &operator,
            &operator,
            &workspace.rhs,
            controls.relative_tolerance,
            1,
            &mut || !operator.poll(),
        );
        report.linear_iterations += linear.iters - before;
        if let Some(error) = error.borrow_mut().take() {
            return Err(error.into());
        }
        if outcome.progress != SolverRunProgress::Complete {
            return Err(ImplicitConductionError::Invalid(
                "interrupted conduction Krylov callback",
            ));
        }
        if outcome.report.converged {
            workspace.direction.copy_from_slice(&linear.x);
            return Ok(());
        }
        if outcome.report.diagnosis == Some(fs_solver::StallDiagnosis::Breakdown) {
            break;
        }
    }
    Err(ImplicitConductionError::KrylovNotConverged {
        relative_residual: linear.rel_residual(),
    })
}

#[cfg(test)]
#[path = "graph/tests.rs"]
mod tests;
