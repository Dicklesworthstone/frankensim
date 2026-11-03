//! Coupled damped Newton solves for source-backed steady contacts.
//!
//! The physical residual and bulk tangent remain owned by the original
//! conduction assemblers. The contact owner supplies only its additional
//! dK_contact/dT * T feedback. The numerical iteration, inexact FGMRES forcing,
//! Armijo damping and iteration report remain owned by fs-solver.

use std::cell::{Cell, RefCell};

use fs_exec::Cx;
use fs_solver::{
    Globalization, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    NewtonStallDiagnosis, NonlinearProblem, SolverRunProgress, norm2,
};
use fs_sparse::precond::Precond;
use fs_sparse::{Coo, Csr};

use crate::assemble::{
    AssembledSystem, assemble_jacobian_with_optional_interfaces,
    assemble_operator_scaled_with_interfaces, full_residual, reduce_matrix_and_lift,
};
use crate::solve::{energy_balance, spd_preconditioner};
use crate::{
    ConductionError, ConductionProblem, DofMap, EnergyBalance, RobinFlux,
};

use super::super::{TemperatureContactEvaluation, TemperatureContacts, invalid, poll, positive};

/// Independent physical acceptance gates and bounded Newton work.
///
/// The numerical absolute tolerance may be tighter than the watt residual
/// gate, but never looser. Relative stopping and trust-region globalization
/// are deliberately not admitted by this damped-Newton entry point.
#[derive(Debug, Clone, Copy)]
pub struct SteadyContactNewtonConfig {
    /// Original numerical owner's forcing, Krylov and Armijo controls.
    pub newton: NewtonKrylovConfig,
    /// Maximum attempted outer iterations; positive.
    pub max_iterations: usize,
    /// Maximum physical residual callbacks, including the initial field and
    /// rejected line-search trials. Tangent assemblies are counted separately.
    pub max_residual_evaluations: usize,
    /// Original contact tangent's matching-pair/retained-triangle allowance.
    /// The contact set's existing independent quadrature allowance also applies.
    pub max_face_pairs: usize,
    /// Fixed absolute Euclidean norm limit on all free residual rows, W.
    pub residual_tolerance_w: f64,
    /// Fixed absolute independent global heat-closure limit, W.
    pub energy_tolerance_w: f64,
}

impl Default for SteadyContactNewtonConfig {
    fn default() -> Self {
        Self {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1.0e-8,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            max_iterations: 32,
            max_residual_evaluations: 4096,
            max_face_pairs: 1_000_000,
            residual_tolerance_w: 1.0e-8,
            energy_tolerance_w: 1.0e-8,
        }
    }
}

impl SteadyContactNewtonConfig {
    fn admit(self, dimension: usize) -> Result<(), ConductionError> {
        self.newton.validate().map_err(|error| {
            invalid("<steady-contact-newton>", format!("Newton configuration refused: {error}"))
        })?;
        if self.max_iterations == 0 || self.max_residual_evaluations == 0
            || self.max_face_pairs == 0
            || !(self.residual_tolerance_w.is_finite() && self.residual_tolerance_w > 0.0)
            || !(self.energy_tolerance_w.is_finite() && self.energy_tolerance_w > 0.0)
            || self.newton.relative_tolerance != 0.0
            || self.newton.absolute_tolerance > self.residual_tolerance_w
            || !matches!(self.newton.globalization, Globalization::LineSearch(_))
        {
            return Err(invalid("<steady-contact-newton>",
                "positive work and absolute watt budgets, line-search globalization, zero relative tolerance and a numerical target no looser than the physical target are required"));
        }
        self.newton.linear_restart.checked_mul(2)
            .and_then(|work| work.checked_add(2))
            .and_then(|work| work.checked_mul(self.newton.max_linear_cycles))
            .and_then(|work| work.checked_mul(self.max_iterations))
            .and_then(|work| work.checked_add(self.max_residual_evaluations))
            .ok_or_else(|| invalid("<steady-contact-newton>", "combined Newton work allowance overflows"))?;
        // Check the original Krylov workspace arithmetic before its allocation.
        let restart = self.newton.linear_restart;
        restart.checked_add(1).and_then(|rows| {
            let basis = rows.checked_mul(dimension)?.checked_mul(2)?;
            let hessenberg = rows.checked_mul(restart)?;
            let small = rows.checked_mul(3)?;
            basis.checked_add(hessenberg)?.checked_add(small)?
                .checked_mul(core::mem::size_of::<f64>())
        }).filter(|&bytes| bytes <= isize::MAX as usize)
            .ok_or_else(|| invalid("<steady-contact-newton>", "Krylov workspace allowance overflows"))?;
        Ok(())
    }
}

/// An accepted physical field and the numerical owner's unmodified report.
///
/// This is not a frozen-contact ConductionSolution: there was no replacement
/// final PDE solve, and no linear-solve report is invented for this result.
#[derive(Debug)]
pub struct ContactSteadyNewtonSolution {
    /// Complete positive Kelvin field, with prescribed entries exact.
    pub temperature: Vec<f64>,
    /// Source receipts queried at this very field, preserving original pins.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Original global vertex ids and physical residual values, W.
    pub free_residual_w: Vec<(usize, f64)>,
    /// Euclidean norm of the complete free-row residual, W.
    pub residual_norm_w: f64,
    /// Source, boundary, reservoir and global closure accounting.
    pub energy: EnergyBalance,
    /// Signed outward Robin transfers at the accepted field.
    pub robin_fluxes: Vec<RobinFlux>,
    /// Original fixed vertex ids and imposed reservoir input reactions, W.
    pub dirichlet_reactions_w: Vec<(usize, f64)>,
    /// Original fs-solver report, including forcing and accepted damping.
    pub newton: NewtonReport,
    /// Actual attempted residual callbacks, including rejected source trials.
    pub residual_evaluations: usize,
    /// Complete tangents assembled; reused for every Krylov action at that field.
    pub jacobian_assemblies: usize,
}

struct PhysicalEvaluation {
    free_temperature: Vec<f64>,
    temperature: Vec<f64>,
    contacts: TemperatureContactEvaluation,
    system: AssembledSystem,
    residual: Vec<f64>,
}

struct Tangent {
    free_temperature: Vec<f64>,
    matrix: Csr,
    preconditioner: Box<dyn Precond>,
}

/// All mutation is private scratch; no source card or caller field is mutated.
/// fs-solver callbacks have no Result channel, so failures are retained here.
/// Only admitted trial-domain misses become NaNs that the original Armijo
/// driver rejects. Fatal failures also trip its transactional control poll.
struct PhysicalProblem<'a, 'm, 'c, 'cx> {
    cx: &'a Cx<'cx>,
    contacts: &'a TemperatureContacts<'m, 'c>,
    problem: ConductionProblem<'a>,
    dofs: &'a DofMap,
    config: SteadyContactNewtonConfig,
    physical: RefCell<Option<PhysicalEvaluation>>,
    tangent: RefCell<Option<Tangent>>,
    fatal: RefCell<Option<ConductionError>>,
    trial_refusal: RefCell<Option<ConductionError>>,
    residual_evaluations: Cell<usize>,
    jacobian_assemblies: Cell<usize>,
}

fn same_field(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}

fn finite(value: f64, field: &'static str) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(ConductionError::NonFinite { field, bits: value.to_bits() }) }
}

enum EvaluationFailure {
    Domain(ConductionError),
    Fatal(ConductionError),
}

impl From<ConductionError> for EvaluationFailure {
    fn from(error: ConductionError) -> Self {
        if matches!(error, ConductionError::OutsideTemperatureSpan { .. }) {
            Self::Domain(error)
        } else {
            Self::Fatal(error)
        }
    }
}

impl PhysicalProblem<'_, '_, '_, '_> {
    fn fail(&self, error: ConductionError) {
        let mut fatal = self.fatal.borrow_mut();
        if fatal.is_none() { *fatal = Some(error); }
    }

    fn stopped(&self) -> bool {
        if self.fatal.borrow().is_some() { return true; }
        if let Err(error) = poll(self.cx, self.residual_evaluations.get()) {
            self.fail(error);
            return true;
        }
        false
    }

    fn record_residual_failure(&self, failure: EvaluationFailure) {
        match failure {
            EvaluationFailure::Domain(error) => *self.trial_refusal.borrow_mut() = Some(error),
            EvaluationFailure::Fatal(error) => self.fail(error),
        }
    }

    fn evaluate(&self, free_temperature: &[f64]) -> Result<PhysicalEvaluation, EvaluationFailure> {
        let temperature = self.dofs.scatter(free_temperature);
        let mut source_domain = false;
        let contacts = self.contacts.evaluate_with_query_refusal(
            self.cx, self.problem.boundary, &temperature, &mut |error| {
                source_domain = matches!(error,
                    // The original query owner checks the shared domain set
                    // before selecting the pinned claim. An empty set is
                    // therefore also a source-domain miss for a pinned query.
                    fs_matdb::MatDbError::NoClaimInDomain { .. }
                    | fs_matdb::MatDbError::PinnedClaimOutOfDomain { .. }
                    | fs_matdb::MatDbError::OutsideKnotSpan { .. }
                    | fs_matdb::MatDbError::OutsideSurfaceSupport { .. });
            },
        ).map_err(|error| if source_domain {
            EvaluationFailure::Domain(error)
        } else {
            EvaluationFailure::Fatal(error)
        })?;
        let system = assemble_operator_scaled_with_interfaces(
            self.cx, self.problem.mesh, self.problem.boundary, self.problem.material,
            self.problem.source, &temperature, None, Some(contacts.frozen_interfaces()),
            self.problem.element_materials,
        )?;
        let residual = full_residual(&system, &temperature);
        for (i, &value) in residual.iter().enumerate() {
            if i % 256 == 0 { poll(self.cx, i)?; }
            finite(value, "steady contact physical residual")?;
        }
        Ok(PhysicalEvaluation {
            free_temperature: free_temperature.to_vec(), temperature, contacts, system, residual,
        })
    }

    fn prepare_tangent(&self, free_temperature: &[f64]) -> Result<(), ConductionError> {
        if self.tangent.borrow().as_ref()
            .is_some_and(|tangent| same_field(&tangent.free_temperature, free_temperature))
        {
            return Ok(());
        }
        poll(self.cx, self.jacobian_assemblies.get())?;
        let count = self.jacobian_assemblies.get().checked_add(1)
            .filter(|&count| count <= self.config.max_iterations)
            .ok_or_else(|| invalid("<steady-contact-newton>", "complete tangent assembly allowance exhausted"))?;
        self.jacobian_assemblies.set(count);
        // A previous accepted step may have rejected trials before acceptance.
        // Preserve only refusals from the current outer attempt on failure.
        self.trial_refusal.borrow_mut().take();
        let physical = self.physical.borrow();
        let physical = physical.as_ref().filter(|physical|
            same_field(&physical.free_temperature, free_temperature))
            .ok_or_else(|| invalid("<steady-contact-newton>", "tangent field has no matching physical residual"))?;
        let contact = self.contacts.linearize_temperature_for_newton(
            self.cx, self.problem.boundary, &physical.temperature, self.config.max_face_pairs,
        )?;
        let bulk_and_frozen = assemble_jacobian_with_optional_interfaces(
            self.cx, self.problem.mesh, self.problem.boundary, self.problem.material,
            &physical.temperature, Some(contact.physical_contacts().frozen_interfaces()),
            self.problem.element_materials,
        )?;
        let n = self.problem.mesh.vertex_count();
        let mut complete = Coo::new(n, n);
        // bulk_and_frozen ALREADY includes K_contact. feedback contains only
        // dK_contact/dT * T, so the physical contact block is counted once.
        for row in 0..n {
            if row % 256 == 0 { poll(self.cx, row)?; }
            for matrix in [&bulk_and_frozen, contact.feedback()] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    complete.push(row, column, finite(value, "steady contact tangent entry")?);
                }
            }
        }
        let complete = complete.assemble();
        for row in 0..n {
            if row % 256 == 0 { poll(self.cx, row)?; }
            for &value in complete.row(row).1 {
                finite(value, "assembled steady contact tangent entry")?;
            }
        }
        // Fixed perturbations are zero. The primal Dirichlet lift is not a
        // tangent load. Precondition with the original frozen physical owner.
        let (matrix, _) = reduce_matrix_and_lift(&complete, self.dofs);
        let (physical_matrix, _) = reduce_matrix_and_lift(&physical.system.operator, self.dofs);
        let preconditioner: Box<dyn Precond> = Box::new(spd_preconditioner(&physical_matrix));
        poll(self.cx, count)?;
        *self.tangent.borrow_mut() = Some(Tangent {
            free_temperature: free_temperature.to_vec(), matrix, preconditioner,
        });
        Ok(())
    }
}

impl NonlinearProblem for PhysicalProblem<'_, '_, '_, '_> {
    fn dimension(&self) -> usize { self.dofs.n() }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.stopped() { return; }
        if x.len() != self.dimension() || output.len() != self.dimension() {
            self.fail(invalid("<steady-contact-newton>", "residual callback layout differs from the fixed free map"));
            return;
        }
        let Some(count) = self.residual_evaluations.get().checked_add(1)
            .filter(|&count| count <= self.config.max_residual_evaluations) else {
                self.fail(invalid("<steady-contact-newton>", "physical residual evaluation allowance exhausted"));
                return;
            };
        self.residual_evaluations.set(count);
        for (i, &temperature) in x.iter().enumerate() {
            if i % 256 == 0
                && let Err(error) = poll(self.cx, i)
            {
                self.fail(error);
                return;
            }
            if let Err(error) = finite(temperature, "steady contact trial temperature") {
                self.fail(error);
                return;
            }
            if temperature <= 0.0 {
                *self.trial_refusal.borrow_mut() = Some(invalid("<steady-contact-newton>",
                    "trial absolute temperature is outside the positive Kelvin domain"));
                return;
            }
        }
        match self.evaluate(x) {
            Ok(physical) => {
                for (value, &vertex) in output.iter_mut().zip(self.dofs.free()) {
                    *value = physical.residual[vertex];
                }
                if let Err(error) = finite(norm2(output), "steady contact residual norm") {
                    output.fill(f64::NAN);
                    self.fail(error);
                    return;
                }
                *self.physical.borrow_mut() = Some(physical);
            }
            Err(error) => self.record_residual_failure(error),
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.stopped() { return; }
        if x.len() != self.dimension() || direction.len() != self.dimension()
            || output.len() != self.dimension()
        {
            self.fail(invalid("<steady-contact-newton>", "Jacobian callback layout differs from the fixed free map"));
            return;
        }
        if let Err(error) = self.prepare_tangent(x) {
            self.fail(error);
            return;
        }
        if let Some(tangent) = self.tangent.borrow().as_ref() {
            tangent.matrix.spmv(direction, output);
        }
        for &value in output.iter() {
            if let Err(error) = finite(value, "steady contact Jacobian action") {
                self.fail(error);
                return;
            }
        }
    }

    fn preconditioner_apply(
        &self, x: &[f64], _outer_iteration: usize, _inner_iteration: usize,
        residual: &[f64], output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        if self.stopped() { return; }
        if x.len() != self.dimension() || residual.len() != self.dimension()
            || output.len() != self.dimension()
        {
            self.fail(invalid("<steady-contact-newton>", "preconditioner layout differs from the fixed free map"));
            return;
        }
        if let Err(error) = self.prepare_tangent(x) {
            self.fail(error);
            return;
        }
        if let Some(tangent) = self.tangent.borrow().as_ref() {
            tangent.preconditioner.apply(residual, output);
        }
        for &value in output.iter() {
            if let Err(error) = finite(value, "steady contact preconditioner action") {
                self.fail(error);
                return;
            }
        }
    }
}

impl TemperatureContacts<'_, '_> {
    /// Solve the complete steady contact residual by damped Newton--Krylov.
    ///
    /// Bulk k(T), heterogeneous assignments, matching face means and retained
    /// nonmatching quadrature all use their original owners. The complete
    /// contact derivative is assembled once per Newton field and cached for
    /// all matrix actions. The original frozen physical preconditioner is
    /// reused; no dense Jacobian, finite difference or overlap search is added.
    ///
    /// The initial field is explicit positive Kelvin. Its fixed entries are
    /// replaced by exact prescribed temperatures before the first query.
    /// Source-domain misses shrink the original Armijo step without switching
    /// pinned claims or extrapolating. Other source refusals, arithmetic and
    /// cancellation abort. A rejected final line search retains its original
    /// source-domain refusal when one occurred.
    ///
    /// Publication requires the original numerical convergence report AND the
    /// unchanged absolute physical residual and heat-closure gates. The cached
    /// physical evaluation must match the accepted field bit for bit; it is
    /// not a second PDE solve. Every refusal returns no field and mutates no
    /// caller input. Continuous source curves use the original primal tangent
    /// at knots and inward source-domain endpoints; surface laws retain their
    /// unique-gradient requirement. This forward-solve tangent does not grant
    /// a two-sided endpoint derivative or adjoint, or establish convergence or
    /// uniqueness for arbitrary source laws.
    #[allow(clippy::too_many_lines)]
    pub fn solve_steady_newton(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial_temperature: &[f64],
        config: SteadyContactNewtonConfig,
    ) -> Result<ContactSteadyNewtonSolution, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        let n = problem.mesh.vertex_count();
        if initial_temperature.len() != n {
            return Err(invalid("<steady-contact-newton>", "initial field and contact mesh differ"));
        }
        for (i, &temperature) in initial_temperature.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            positive(temperature)?;
        }
        problem.source.validate("volumetric source", n)?;
        if let Some(assigned) = problem.element_materials { assigned.validate_for(problem.mesh)?; }
        let dofs = DofMap::new(problem.boundary, n)?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        config.admit(dofs.n())?;
        let initial = dofs.gather(initial_temperature);
        for (i, &vertex) in dofs.fixed().iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            positive(dofs.prescribed()[vertex])?;
        }
        let adapter = PhysicalProblem {
            cx, contacts: self, problem, dofs: &dofs, config,
            physical: RefCell::new(None), tangent: RefCell::new(None),
            fatal: RefCell::new(None), trial_refusal: RefCell::new(None),
            residual_evaluations: Cell::new(0), jacobian_assemblies: Cell::new(0),
        };
        let initialized = NewtonKrylovState::new(&adapter, initial, config.newton);
        if let Some(error) = adapter.fatal.borrow_mut().take() { return Err(error); }
        let mut state = match initialized {
            Ok(state) => state,
            Err(error) => {
                if let Some(refusal) = adapter.trial_refusal.borrow_mut().take() {
                    return Err(refusal);
                }
                return Err(invalid("<steady-contact-newton>", format!("Newton initialization refused: {error}")));
            }
        };
        let mut stopped = || adapter.stopped();
        let run = state.run_controlled(&adapter, config.max_iterations, &mut stopped);
        if let Some(error) = adapter.fatal.borrow_mut().take() { return Err(error); }
        if !matches!(run.progress, SolverRunProgress::Complete) {
            return Err(invalid("<steady-contact-newton>",
                format!("controlled Newton run refused: {:?}", run.progress)));
        }
        let report = run.report;
        if !report.converged {
            if matches!(report.diagnosis, Some(NewtonStallDiagnosis::LineSearchRejected | NewtonStallDiagnosis::NonFinite))
                && let Some(error) = adapter.trial_refusal.borrow_mut().take()
            {
                return Err(error);
            }
            return Err(invalid("<steady-contact-newton>", format!(
                "coupled Newton refused after {} iterations: {:?}; physical residual {} W, target {} W",
                report.iterations, report.diagnosis, report.residual_norm, config.residual_tolerance_w,
            )));
        }
        poll(cx, report.iterations)?;
        let physical = adapter.physical.borrow_mut().take()
            .filter(|physical| same_field(&physical.free_temperature, &state.x))
            .ok_or_else(|| invalid("<steady-contact-newton>", "accepted Newton field lost its matching physical evaluation"))?;
        let free_residual_w: Vec<_> = dofs.free().iter()
            .map(|&vertex| (vertex, physical.residual[vertex])).collect();
        let values: Vec<_> = free_residual_w.iter().map(|&(_, value)| value).collect();
        let residual_norm_w = finite(norm2(&values), "accepted steady contact residual norm")?;
        let (energy, robin_fluxes) = energy_balance(
            problem.mesh, problem.boundary, problem.source,
            &physical.system, &dofs, &physical.temperature,
        );
        for value in [energy.source_w, energy.neumann_out_w, energy.robin_out_w,
            energy.dirichlet_in_w, energy.closure_w, energy.scale_w]
        {
            finite(value, "accepted steady contact heat balance")?;
        }
        if residual_norm_w > config.residual_tolerance_w
            || energy.closure_w.abs() > config.energy_tolerance_w
        {
            return Err(invalid("<steady-contact-newton>", format!(
                "accepted numerical state fails physical gates: residual {residual_norm_w} W (limit {} W), heat closure {} W (limit {} W)",
                config.residual_tolerance_w, energy.closure_w, config.energy_tolerance_w,
            )));
        }
        let dirichlet_reactions_w = dofs.fixed().iter()
            .map(|&vertex| (vertex, physical.residual[vertex])).collect();
        poll(cx, report.iterations)?;
        Ok(ContactSteadyNewtonSolution {
            temperature: physical.temperature, physical_contacts: physical.contacts,
            free_residual_w, residual_norm_w, energy, robin_fluxes, dirichlet_reactions_w,
            newton: report, residual_evaluations: adapter.residual_evaluations.get(),
            jacobian_assemblies: adapter.jacobian_assemblies.get(),
        })
    }
}
