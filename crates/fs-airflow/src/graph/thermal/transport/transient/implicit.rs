//! Sparse frozen-flow enthalpy transport on the original backward-Euler owner.
//!
//! The primary enthalpy equations are linear for prescribed mass flow even
//! when the original temperature/enthalpy law is nonlinear. Native sparse
//! assembly and cancellable FGMRES replace the small dense feedback solve;
//! original caloric inversion and local/global energy acceptance still decide
//! whether a proposed thermal state can be published.

use fs_exec::Cx;
use fs_solver::{CsrOp, FgmresState, SolveReport};
use fs_sparse::precond::Precond;
use fs_sparse::{Coo, Csr};
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::backward_euler::BackwardEulerStep;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};

use super::{GraphEnthalpyTransport, TransientTransportError, thermal_checkpoint};

pub const MAX_SPARSE_ADVECTION_CELLS: usize = 65_536;
pub const MAX_SPARSE_ADVECTION_ENTRIES: usize = 1_048_576;

/// Explicit spatial, Krylov and physical work limits. Linear convergence is
/// separate from the original per-cell and total joule acceptance thresholds.
#[derive(Clone, Copy, Debug)]
pub struct SparseAdvectionControls {
    pub maximum_cells: usize,
    /// Staged coefficients, before deterministic duplicate accumulation.
    pub maximum_entries: usize,
    /// Native FGMRES basis length, bounded to 64.
    pub restart: usize,
    pub maximum_cycles: usize,
    /// Bounded defect corrections after rounding to representable enthalpy.
    pub maximum_corrections: usize,
    pub relative_linear_tolerance: f64,
    pub cell_tolerance_j: f64,
    pub total_tolerance_j: f64,
}

#[derive(Debug)]
pub struct SparseAdvectionReport {
    pub thermal: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    pub matrix_entries: usize,
    /// Actual native true-Euclidean residual reports for accepted corrections.
    pub linear: Vec<SolveReport>,
}

#[derive(Debug)]
pub enum SparseAdvectionError {
    Invalid(&'static str),
    Transport(TransientTransportError),
    Thermal(DepositionError),
    Linear(SolveReport),
}
impl From<TransientTransportError> for SparseAdvectionError {
    fn from(error: TransientTransportError) -> Self {
        Self::Transport(error)
    }
}
impl From<DepositionError> for SparseAdvectionError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for SparseAdvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sparse implicit airflow transport: {self:?}")
    }
}
impl std::error::Error for SparseAdvectionError {}

struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        for ((out, value), inverse) in output.iter_mut().zip(residual).zip(&self.0) {
            *out = value * inverse;
        }
    }
}

impl GraphEnthalpyTransport {
    /// Advance prescribed-flow storage with sparse native linear work.
    ///
    /// All active cells must belong to this transport binding. Native fixed
    /// sources participate; active conduction, ambient/reservoir exchange and
    /// state-dependent flow require the general coupled endpoint interface.
    /// This restriction preserves an exact enthalpy-space linear operator,
    /// without freezing nonlinear temperature-dependent heat exchange.
    ///
    /// The same native upwind and backward-Euler equations as the small dense
    /// path are solved. Sparse storage is O(cells + connections); the native
    /// Krylov workspace is bounded by the declared cells and restart length.
    /// The matrix includes the derivative of the original component enthalpy
    /// reference, including admitted mass-balance roundoff.
    ///
    /// No thermal clock changes on cancellation, source refusal or exhaustion.
    /// Implicit stability does not certify residence-time resolution: refine
    /// both timestep and the original fluid control volumes.
    pub fn advance_advection_sparse_implicit<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: SparseAdvectionControls,
    ) -> Result<(DepositionSimulation<'m>, SparseAdvectionReport), SparseAdvectionError> {
        self.require_upwind()?;
        self.validate_storage(cx, thermal)?;
        let n = thermal.cells().len();
        if !(1..=MAX_SPARSE_ADVECTION_CELLS).contains(&controls.maximum_cells)
            || n == 0
            || n != self.storage.len()
            || n > controls.maximum_cells
            || !(1..=MAX_SPARSE_ADVECTION_ENTRIES).contains(&controls.maximum_entries)
            || !(1..=64).contains(&controls.restart)
            || !(1..=1024).contains(&controls.maximum_cycles)
            || !(1..=4).contains(&controls.maximum_corrections)
            || !controls.relative_linear_tolerance.is_finite()
            || controls.relative_linear_tolerance <= 0.0
            || controls.relative_linear_tolerance >= 1.0
            || !nonnegative(controls.cell_tolerance_j)
            || !nonnegative(controls.total_tolerance_j)
            || thermal
                .links()
                .values()
                .any(|link| link.conductance_w_k != 0.0)
            || thermal
                .cells()
                .values()
                .any(|cell| cell.boundary.htc_w_m2_k != 0.0 || cell.boundary.reservoir.is_some())
        {
            return Err(SparseAdvectionError::Invalid(
                "bounded mapped cells, frozen advection/fixed sources and explicit solve limits required",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
        let ids: std::collections::BTreeMap<_, _> = step
            .cell_ids()
            .iter()
            .enumerate()
            .map(|(i, &id)| (id, i))
            .collect();
        let mass: Vec<_> = step
            .cell_ids()
            .iter()
            .map(|id| thermal.cells()[id].mass_kg)
            .collect();
        let active = self
            .advection
            .flows()
            .iter()
            .filter(|flow| flow.kg_s != 0.0)
            .count();
        let outlets = self
            .advection
            .boundaries()
            .iter()
            .filter(|b| matches!(b, MassBoundary::Outlet { .. }))
            .count();
        let staged = n
            .checked_mul(2)
            .and_then(|v| active.checked_mul(2).and_then(|a| v.checked_add(a)))
            .and_then(|v| v.checked_add(outlets))
            .ok_or(SparseAdvectionError::Invalid(
                "sparse coefficient count overflow",
            ))?;
        if staged > controls.maximum_entries {
            return Err(SparseAdvectionError::Invalid(
                "sparse coefficient budget exceeded",
            ));
        }
        let dt = target_s - thermal.time_s();
        let mut matrix = Coo::new(n, n);
        let mut balance = vec![0.0; n];
        let mut components: Vec<_> = (0..n).collect();
        for i in 0..n {
            thermal_checkpoint(cx)?;
            matrix.push(i, i, 1.0);
        }
        for flow in self
            .advection
            .flows()
            .iter()
            .filter(|flow| flow.kg_s != 0.0)
        {
            thermal_checkpoint(cx)?;
            let (up, down) = if flow.kg_s > 0.0 {
                (ids[&flow.a], ids[&flow.b])
            } else {
                (ids[&flow.b], ids[&flow.a])
            };
            let rate = flow.kg_s.abs();
            let coefficient = finite(dt * (rate / mass[up]))?;
            matrix.push(up, up, coefficient);
            matrix.push(down, up, -coefficient);
            balance[up] = finite(balance[up] - rate)?;
            balance[down] = finite(balance[down] + rate)?;
            let a = component_root(&mut components, up);
            let b = component_root(&mut components, down);
            components[a.max(b)] = a.min(b);
        }
        for boundary in self.advection.boundaries() {
            thermal_checkpoint(cx)?;
            match *boundary {
                MassBoundary::Inlet { cell, kg_s, .. } => {
                    let i = ids[&cell];
                    balance[i] = finite(balance[i] + kg_s)?;
                }
                MassBoundary::Outlet { cell, kg_s, .. } => {
                    let i = ids[&cell];
                    balance[i] = finite(balance[i] - kg_s)?;
                    matrix.push(i, i, finite(dt * (kg_s / mass[i]))?);
                }
            }
        }
        for i in 0..n {
            thermal_checkpoint(cx)?;
            let reference = component_root(&mut components, i);
            matrix.push(i, reference, finite(dt * (balance[i] / mass[reference]))?);
        }
        thermal_checkpoint(cx)?;
        let matrix: Csr = matrix.assemble();
        thermal_checkpoint(cx)?;
        let mut inverses = Vec::with_capacity(n);
        for i in 0..n {
            thermal_checkpoint(cx)?;
            if matrix.row(i).1.iter().any(|v| !v.is_finite()) {
                return Err(SparseAdvectionError::Invalid(
                    "assembled enthalpy operator is nonfinite",
                ));
            }
            let diagonal = matrix.get(i, i);
            let inverse = 1.0 / diagonal;
            if !diagonal.is_finite() || diagonal <= 0.0 || !inverse.is_finite() || inverse <= 0.0 {
                return Err(SparseAdvectionError::Invalid(
                    "positive representable transport diagonal required",
                ));
            }
            inverses.push(inverse);
        }
        let entries = matrix.nnz();
        let operator = CsrOp::general(matrix);
        let preconditioner = Diagonal(inverses);
        let mut energies = vec![0.0; n];
        let mut reports = Vec::new();
        for correction in 0..=controls.maximum_corrections {
            thermal_checkpoint(cx)?;
            let trial = step.evaluate_with_transport(
                &energies,
                |candidate| {
                    self.advection
                        .endpoint_transport(candidate, || thermal_checkpoint(cx))
                },
                || thermal_checkpoint(cx),
            )?;
            let maximum_cell_residual_j = trial.maximum_cell_residual_j();
            if correction == controls.maximum_corrections
                || (maximum_cell_residual_j <= controls.cell_tolerance_j
                    && trial.audit().residual_j.abs() <= controls.total_tolerance_j)
            {
                let (next, (), audit) = trial.accept(
                    controls.cell_tolerance_j,
                    controls.total_tolerance_j,
                    || thermal_checkpoint(cx),
                )?;
                thermal_checkpoint(cx)?;
                return Ok((
                    next,
                    SparseAdvectionReport {
                        thermal: audit,
                        maximum_cell_residual_j,
                        matrix_entries: entries,
                        linear: reports,
                    },
                ));
            }
            let scale = trial
                .residual_j()
                .iter()
                .map(|r| r.abs())
                .fold(0.0f64, f64::max);
            if scale == 0.0 {
                // A global native refusal cannot be repaired by inventing an
                // external source or changing an already-zero local residual.
                let _ = trial.accept(
                    controls.cell_tolerance_j,
                    controls.total_tolerance_j,
                    || thermal_checkpoint(cx),
                )?;
                return Err(SparseAdvectionError::Invalid(
                    "zero local residual with unresolved global balance",
                ));
            }
            let rhs: Vec<_> = trial.residual_j().iter().map(|r| -r / scale).collect();
            let mut solver = FgmresState::new(&rhs, controls.restart.min(n));
            let run = solver.run_cancellable(
                &operator,
                &preconditioner,
                &rhs,
                controls.relative_linear_tolerance,
                controls.maximum_cycles,
                cx,
            );
            thermal_checkpoint(cx)?;
            if !run.report.converged_euclidean() {
                return Err(SparseAdvectionError::Linear(run.report));
            }
            for (energy, delta) in energies.iter_mut().zip(&solver.x) {
                *energy = finite(*energy + finite(*delta * scale)?)?;
            }
            energies = step.resolvable_increments_j(&energies);
            reports.push(run.report);
        }
        unreachable!("bounded correction loop returns through original thermal acceptance")
    }
}

// Compress paths while preserving the native lowest-CellId representative.
fn component_root(parents: &mut [usize], mut index: usize) -> usize {
    let mut root = index;
    while parents[root] != root {
        root = parents[root];
    }
    while parents[index] != index {
        let next = parents[index];
        parents[index] = root;
        index = next;
    }
    root
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, SparseAdvectionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SparseAdvectionError::Invalid(
            "unrepresentable sparse transport arithmetic",
        ))
    }
}

#[cfg(test)]
mod tests;
