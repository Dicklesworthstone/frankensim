//! Sparse nonlinear heat transport on the original fixed thermal graph.
use super::*;
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolverRunProgress};
use fs_sparse::Csr;
use fs_sparse::precond::{Ilu0, IluControlledError, ilu0_controlled};
use std::cell::RefCell;

#[path = "graph/forest.rs"]
mod forest;

/// Bounded linear work inside the original Newton/time-step controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionLinearControls {
    /// Arnoldi dimension, 1..=128. Storage includes cells * restart, links
    /// and any explicitly bounded coupled-source dependencies.
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
/// as `ImplicitConductionChain`. Bounded FGMRES replaces the chain's O(n)
/// direct solve. Acyclic Jacobians use exact sparse forest preconditioning;
/// cyclic graphs use zero-fill ILU. Unsafe factors retain diagonal scaling.
pub struct ImplicitConductionGraph<'a> {
    inner: ImplicitConductionChain<'a>,
}
impl<'a> ImplicitConductionGraph<'a> {
    /// Solve sparse nonlocal heat-generation feedback with graph conduction.
    /// See the chain attachment for immutable original-caloric bindings,
    /// sparse derivative units, source ownership and physical-time semantics.
    /// Non-neighbor source dependencies use this graph's existing Jv/J^T and
    /// bounded Krylov solve; no alternate thermal integration is introduced.
    pub fn with_coupled_heat_source(
        self,
        thermal: &DepositionSimulation<'_>,
        cells: &[CellId],
        law: &'a dyn super::coupled_source::ImplicitCoupledHeatSourceLaw,
        maximum_derivatives: usize,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_coupled_heat_source(
                thermal, cells, law, maximum_derivatives, checkpoint,
            )?,
        })
    }

    /// Solve original local heat-generation feedback with graph conduction.
    /// The chain attachment documents caloric identity, signed-power units and
    /// additive source ownership. Sparse topology and linear budgets are kept.
    pub fn with_heat_sources(
        self,
        thermal: &DepositionSimulation<'_>,
        sources: &[(CellId, &'a dyn super::source::ImplicitHeatSourceLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_heat_sources(thermal, sources, checkpoint)?,
        })
    }

    /// Resolve two actual prescribed baths through one selected physical face.
    /// The original rate/entropy owner retains each signed bath contribution
    /// through backward Euler and native SDIRK stages. See the chain attachment
    /// for caloric identity, geometry and exclusive channel ownership checks.
    pub fn with_boundary_exchanges(
        self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a dyn super::boundary::ImplicitBoundaryExchangeLaw)],
        ambient: &[(CellId, &'a dyn super::boundary::ImplicitBoundaryExchangeLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_boundary_exchanges(thermal, reservoirs, ambient, checkpoint)?,
        })
    }

    /// Bind original endpoint upwind enthalpy to this same nonlinear graph.
    /// Directed flows may join any active cells. The analytic nonsymmetric
    /// tangent uses the existing sparse FGMRES path. See the chain method for
    /// storage, source ownership and unsupported-carrier restrictions.
    pub fn with_implicit_advection(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &crate::enthalpy::advection::Advected<()>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_implicit_advection(thermal, flow, checkpoint)?,
        })
    }

    /// Bind continuous prescribed pump/inlet histories to the native graph.
    /// See the chain method for fixed-mass, stage-sampling and zero-flow rules.
    pub fn with_scheduled_implicit_advection(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &crate::enthalpy::advection::Advected<()>,
        schedule: &crate::enthalpy::advection::schedule::AdvectionSchedule,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_scheduled_implicit_advection(thermal, flow, schedule, checkpoint)?,
        })
    }

    /// Schedule an already attached upwind or reconstructed flow.
    pub fn with_implicit_advection_schedule(
        self,
        thermal: &DepositionSimulation<'_>,
        schedule: &crate::enthalpy::advection::schedule::AdvectionSchedule,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_implicit_advection_schedule(thermal, schedule, checkpoint)?,
        })
    }

    /// Classify the already attached bulk or liquid upwind flow on the native entropy
    /// account. See [`ImplicitConductionChain::with_implicit_advection_entropy`]
    /// for the caloric, stage-quadrature and reconstruction restrictions.
    pub fn with_implicit_advection_entropy(
        self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_implicit_advection_entropy(thermal, checkpoint)?,
        })
    }

    /// Solve limited axial enthalpy transport with the original conduction
    /// residual, nonlinear iteration and energy accounts.
    ///
    /// Every trial reconstructs the actual endpoint carried h on the supplied
    /// original geometry. A liquid carrier uses liquid h and its own phase
    /// derivative, preserving stationary-solid latent energy. Its analytic
    /// active-minmod tangent joins the native upwind
    /// tangent in the sparse nonsymmetric solver; the wider stencil requires
    /// this graph path even for a straight channel. At limiter switches Newton
    /// uses a deterministic generalized derivative and the existing line search.
    ///
    /// Empty plans preserve the upwind path. The original reconstruction's
    /// source, topology, fixed-mass, geometry and unsupported-carrier refusals
    /// remain active. Select backward Euler or native SDIRK2 for first- or
    /// second-order time integration, with their existing adaptive controllers.
    /// Limiter switches can reduce observed order; neither a large accepted step
    /// nor a small energy residual certifies transport accuracy or a maximum
    /// principle for arbitrary coupled sources.
    pub fn with_implicit_advection_reconstruction(
        self,
        thermal: &DepositionSimulation<'_>,
        flow: &crate::enthalpy::advection::Advected<()>,
        reconstruction: &crate::enthalpy::advection::reconstruction::LimitedEnthalpyReconstruction,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_implicit_advection_reconstruction(
                thermal,
                flow,
                reconstruction,
                checkpoint,
            )?,
        })
    }

    /// Replace selected reservoir and ambient fluxes in the same original
    /// residual and Jacobian. See the chain method for original caloric,
    /// geometric and boundary-channel ownership requirements.
    pub fn with_boundary_laws(
        self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a dyn super::boundary::ImplicitBoundaryLaw)],
        ambient: &[(CellId, &'a dyn super::boundary::ImplicitBoundaryLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self
                .inner
                .with_boundary_laws(thermal, reservoirs, ambient, checkpoint)?,
        })
    }

    /// Prescribe continuous bath temperatures at each native stage time.
    /// The chain attachment retains original face/caloric ownership, clips
    /// intervals at history knots and accounts heat at the actual bath Kelvin.
    pub fn with_bath_temperature_histories(
        self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a super::boundary::BathTemperatureSchedule)],
        ambient: &[(CellId, &'a super::boundary::BathTemperatureSchedule)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_bath_temperature_histories(
                thermal, reservoirs, ambient, checkpoint,
            )?,
        })
    }

    /// Replace selected bulk-series links with original integrated conductive
    /// fluxes and their full tangents. The same attachment and ownership rules
    /// as [`ImplicitConductionChain::with_link_laws`] apply.
    pub fn with_link_laws(
        self,
        thermal: &DepositionSimulation<'_>,
        links: &[(u64, &'a dyn super::link::ImplicitLinkLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_link_laws(thermal, links, checkpoint)?,
        })
    }

    /// Attach original two-sided contact laws to this same graph and nonlinear
    /// residual. See the chain entry for binding and no-double-counting rules.
    pub fn with_contacts(
        self,
        thermal: &DepositionSimulation<'_>,
        contacts: &[(u64, &'a dyn super::contact::ImplicitContactLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Ok(Self {
            inner: self.inner.with_contacts(thermal, contacts, checkpoint)?,
        })
    }

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

    /// Solve stationary original heat balance with the same sparse tangent,
    /// material bindings and Krylov budgets as native graph evolution.
    pub fn solve_operating_point<'s, 'm>(
        &'s self,
        initial_guess: &DepositionSimulation<'m>,
        controls: super::equilibrium::ThermalOperatingControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<super::equilibrium::ThermalOperatingPoint<'s, 'a, 'm>,
        super::equilibrium::ThermalOperatingError> {
        self.inner.solve_operating_point(initial_guess, controls, workspace, checkpoint)
    }

    /// Prepare the same original heat-rate tangent for continuous thermal
    /// impedance. See the chain method for equilibrium and frozen-clock scope.
    pub fn prepare_thermal_impedance(
        &self,
        thermal: &DepositionSimulation<'_>,
        budget: super::harmonic::ThermalImpedanceBudget,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<super::harmonic::ThermalImpedanceModel, super::harmonic::ThermalImpedanceError> {
        self.inner.prepare_thermal_impedance(thermal, budget, checkpoint)
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

    /// Second-order L-stable convection/conduction on the original graph.
    /// The chain entry point documents the fixed-cell boundary schedules and
    /// original domain/entropy/residual gates.
    /// Attached prescribed advection, including limited reconstruction, uses
    /// each actual stage. Both stages and rejected intervals share this graph's
    /// bounded sparse linear owner.
    pub fn advance_sdirk2_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
        self.inner
            .advance_sdirk2_to(thermal, target_s, controls, workspace, checkpoint)
    }

    /// IMEX sources and conservative transfers through the same graph solve.
    /// Sources and directed transfers are held at each accepted start; the
    /// original flow timestep limit remains active. See the chain's
    /// [`ImplicitConductionChain::advance_observed_to`] for history, retry and
    /// unsupported-hook semantics. No explicit conduction CFL is imposed.
    pub fn advance_observed_to<'m, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'m>,
        observer: &O,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, O, ImplicitConductionReport),
        ImplicitConductionObservedError<O::Error>,
    > {
        self.inner
            .advance_observed_to(thermal, observer, target_s, controls, workspace, checkpoint)
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

    /// Native SDIRK2 with the shared h/T/phase error controller and sparse
    /// linear work budget, including attached prescribed convection. The entire
    /// window must satisfy the chain method's boundary-scheduled, fixed-cell
    /// scope before any event processing.
    pub fn advance_sdirk2_adaptive_to<'m>(
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
            .advance_sdirk2_adaptive_to(thermal, target_s, controls, workspace, checkpoint)
    }
}

// The physical Jacobian has one diagonal and two off-diagonals per link, even
// with unequal masses, nonlinear k or latent plateaus (generally nonsymmetric).
// No dense matrix, symmetry assumption or numerical constitutive difference.
enum SparsePreconditioner {
    Forest(forest::Forest),
    Incomplete(Ilu0),
}

// Canonical assembly for the existing nonsymmetric sparse ILU owner. Map
// insertion keeps duplicate summation in the original contribution order and
// gives bounded checkpoint gaps without an uninterruptible global sort.
fn incomplete_factor(
    diagonal: &[f64],
    entries: &[(usize, usize, f64)],
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Option<Ilu0>, DepositionError> {
    if entries.is_empty() {
        return Ok(None);
    }
    let n = diagonal.len();
    let mut combined = BTreeMap::new();
    for (i, &value) in diagonal.iter().enumerate() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        combined.insert((i, i), value);
    }
    for (i, &(row, column, value)) in entries.iter().enumerate() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        if row >= n || column >= n {
            return Ok(None);
        }
        let entry = combined.entry((row, column)).or_insert(0.0);
        *entry += value;
        if !entry.is_finite() {
            return Ok(None);
        }
    }
    let row_count = n.checked_add(1).ok_or(DepositionError::Budget(
        "conduction sparse factor allocation",
    ))?;
    let mut rows = Vec::new();
    let mut columns = Vec::new();
    let mut values = Vec::new();
    rows.try_reserve_exact(row_count)
        .map_err(|_| DepositionError::Budget("conduction sparse factor allocation"))?;
    rows.resize(row_count, 0);
    columns
        .try_reserve_exact(combined.len())
        .map_err(|_| DepositionError::Budget("conduction sparse factor allocation"))?;
    values
        .try_reserve_exact(combined.len())
        .map_err(|_| DepositionError::Budget("conduction sparse factor allocation"))?;
    for (i, ((row, column), value)) in combined.into_iter().enumerate() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        if row != column && value == 0.0 {
            continue;
        }
        rows[row + 1] += 1;
        columns.push(column);
        values.push(value);
    }
    for i in 1..rows.len() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        rows[i] += rows[i - 1];
    }
    let matrix =
        Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, &mut *checkpoint)?
            .ok_or(DepositionError::Arithmetic)?;
    match ilu0_controlled(&matrix, checkpoint) {
        Ok(factor) => Ok(Some(factor)),
        Err(IluControlledError::Breakdown(_)) => Ok(None),
        Err(IluControlledError::Capacity) => Err(DepositionError::Budget(
            "conduction sparse factor allocation",
        )),
        Err(IluControlledError::Checkpoint(error)) => Err(error),
    }
}

struct GraphOperator<'a, 'c> {
    diagonal: &'a [f64],
    entries: &'a [(usize, usize, f64)],
    preconditioner: Option<SparsePreconditioner>,
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
        let applied = match &self.preconditioner {
            Some(SparsePreconditioner::Forest(forest)) => {
                forest.apply(residual, output, || self.poll())
            }
            Some(SparsePreconditioner::Incomplete(factor)) => {
                match factor.apply_controlled(residual, output, || {
                    if self.poll() { Ok(()) } else { Err(()) }
                }) {
                    Ok(()) => true,
                    Err(IluControlledError::Breakdown(_)) => false,
                    Err(IluControlledError::Capacity) => {
                        *self.error.borrow_mut() = Some(DepositionError::Budget(
                            "conduction sparse factor allocation",
                        ));
                        false
                    }
                    Err(IluControlledError::Checkpoint(())) => false,
                }
            }
            None => false,
        };
        if applied || self.error.borrow().is_some() {
            return;
        }
        // Unsafe preconditioner arithmetic does not reject the original
        // system. Retry this action with the existing diagonal scaling.
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
    if report.linear_iterations >= controls.maximum_iterations {
        return Err(DepositionError::Budget("conduction Krylov iterations").into());
    }
    let entries = workspace
        .graph_entries
        .as_deref()
        .ok_or(ImplicitConductionError::Invalid("missing graph Jacobian"))?;
    let forest =
        forest::Forest::factor(&workspace.diagonal, entries, checkpoint).map_err(|error| {
            match error {
                forest::FactorError::Capacity => {
                    DepositionError::Budget("conduction forest allocation")
                }
                forest::FactorError::Checkpoint(error) => error,
            }
        })?;
    let preconditioner = if let Some(forest) = forest {
        Some(SparsePreconditioner::Forest(forest))
    } else {
        incomplete_factor(&workspace.diagonal, entries, checkpoint)?
            .map(SparsePreconditioner::Incomplete)
    };
    let error = RefCell::new(None);
    let checkpoint: &mut dyn FnMut() -> Result<(), DepositionError> = checkpoint;
    let checkpoint = RefCell::new(checkpoint);
    let operator = GraphOperator {
        diagonal: &workspace.diagonal,
        entries,
        preconditioner,
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
