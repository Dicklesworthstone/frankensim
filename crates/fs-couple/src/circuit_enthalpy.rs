//! Transient Kirchhoff circuits coupled to the shared enthalpy network.
//!
//! Every material resistor has an explicit receiving thermal cell. Conduction,
//! reservoirs, latent storage and local energy residuals stay with fs-thermal;
//! charges, fluxes, KCL/KVL and branch losses stay with MaterialCircuit. No trial
//! advances either accepted history. This is backward Euler for thermal rates
//! with circuit midpoint currents and ENDPOINT material coefficients: first-order
//! coupled accuracy in general, not a second-order or unconditional-stability claim.

/// Material field/current response coupled to the same finite enthalpy solve.
pub mod conductivity;
/// Spatial DC current and local Joule heating on the same finite thermal solve.
pub mod conduction;
/// Data-only caloric/card, resistor and conductive-link bindings.
pub mod material;
/// Two-sided deposition of source-backed electrical contact heat.
pub mod contacts;
/// Reciprocal material Seebeck voltage and Peltier heat on the same network.
pub mod thermoelectric;
/// Joint material electric/magnetic gap motion and finite-body Joule heating.
pub mod motion;
/// Distributed sourced winding loss on the same native thermal residual.
pub mod winding;
/// Reversible sourced C(T)/L(T) feedback on the same finite caloric stores.
pub mod field_storage;
/// Shared original material/volume correspondence for magnetic heat consumers.
pub mod magnetic;
/// Retained spatial induction and Joule heat on the original enthalpy timeline.
pub mod eddy;
/// Piezoelectric motion, electrical load and finite resistor heating.
pub mod piezoelectric;

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;
use crate::material_circuit::{MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitError,
    MaterialCircuitSolve, MaterialCircuitStep, MaterialResistor};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};

/// One complete branch-to-cell heat destination; original graph branch indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResistorCell {
    /// Original circuit branch, not resistor-only order.
    pub branch: usize,
    /// Active receiving thermal control volume.
    pub cell: CellId,
}

/// One electrical history, held source ports and complete heat destinations.
#[derive(Clone, Copy)]
pub struct CircuitEnthalpyInput<'a> {
    /// Immutable accepted charge/flux state; never advanced inside a callback.
    pub circuit: &'a MaterialCircuit,
    /// Held external ports; empty for a source-free finite electrical store.
    pub ports: &'a [f64],
    /// Exactly one destination for every resistor, in any order. Several branches
    /// may feed one cell only when that homogenization is physically justified.
    pub destinations: &'a [ResistorCell],
}

/// One repeatable material evaluation at the actual proposed thermal endpoint.
pub struct CircuitThermalCoefficients<S> {
    /// Complete source-bound resistor states for the existing circuit owner.
    pub resistors: Vec<MaterialResistor>,
    /// Independent signed cell powers [W], ADDITIONAL to authored boundaries.
    /// Do not include circuit dissipation here: it is deposited exactly once.
    pub sources_w: BTreeMap<CellId, f64>,
    /// Explicit replacement conductances [W/K] of existing thermal links.
    pub conductances_w_k: BTreeMap<u64, f64>,
    /// Original material/thermal-interface results, retained without type erasure.
    pub source: S,
}

/// Independent nonlinear, circuit and physical energy/work admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircuitEnthalpySolve {
    /// Maximum active thermal cells, including passive cells outside the circuit.
    pub maximum_cells: usize,
    /// All circuit solves, including Jacobian stencils and final fresh admission.
    pub maximum_circuit_trials: usize,
    /// Existing Newton iteration budget; its inner work is bounded by `newton`.
    pub maximum_newton_iterations: usize,
    /// Maximum-component numerical Jacobian perturbation [J of cell energy].
    pub difference_step_j: f64,
    /// Absolute largest per-cell thermal residual [J].
    pub cell_energy_tolerance_j: f64,
    /// Absolute global thermal residual [J], independent of local cancellation.
    pub thermal_energy_tolerance_j: f64,
    /// Absolute deposited-Joule and combined circuit/thermal energy tolerance [J].
    pub total_energy_tolerance_j: f64,
    /// Existing electrical branch-power and supply/storage admission.
    pub circuit: MaterialCircuitSolve,
    /// Existing Newton--Krylov policy, acting on joule coordinates.
    pub newton: NewtonKrylovConfig,
}

/// Original source errors survive infallible numerical callbacks.
#[derive(Debug)]
pub enum CircuitEnthalpyError<E> {
    /// Incomplete mapping, bad clocks/controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Material or interface callback refusal.
    Source(E),
    /// Existing circuit/storage/branch-work refusal.
    Circuit(MaterialCircuitError),
    /// Existing thermal/phase/topology refusal.
    Thermal(DepositionError),
    /// Existing nonlinear construction refusal.
    Newton(NewtonError),
    /// Existing nonlinear solve exhausted its work or stalled.
    NotConverged(NewtonReport),
    /// Interrupted numerical work; no physical endpoint is returned.
    Interrupted(SolverRunProgress),
    /// Circuit evaluation budget exhausted, including final admission work.
    CircuitBudget,
    /// Cancellation before joint publication.
    Cancelled,
}
impl<E> From<DepositionError> for CircuitEnthalpyError<E> {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled { Self::Cancelled } else { Self::Thermal(e) }
    }
}
impl<E> From<MaterialCircuitError> for CircuitEnthalpyError<E> {
    fn from(e: MaterialCircuitError) -> Self {
        match e { MaterialCircuitError::Cancelled => Self::Cancelled, other => Self::Circuit(other) }
    }
}
impl<E> From<NewtonError> for CircuitEnthalpyError<E> {
    fn from(e: NewtonError) -> Self { Self::Newton(e) }
}
impl<E: core::fmt::Debug> core::fmt::Display for CircuitEnthalpyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "circuit/enthalpy coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CircuitEnthalpyError<E> {}
fn numerical<E>(error: FiniteDifferenceError<CircuitEnthalpyError<E>>) -> CircuitEnthalpyError<E> {
    match error { FiniteDifferenceError::Source(e) => e, FiniteDifferenceError::Invalid(s) => CircuitEnthalpyError::Invalid(s) }
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() { Ok(value) } else { Err(DepositionError::Arithmetic) }
}

/// Independently admitted physical work and the exact endpoint source evaluation.
pub struct CircuitEnthalpyReport<S, D = MaterialCircuitStep> {
    /// Original circuit or reciprocal field-motion step. Its physical storage,
    /// external work and branch heat enter the same thermal account.
    pub electrical: D,
    /// Original thermal account; its external heat INCLUDES deposited Joule heat.
    pub thermal: EnergyAudit,
    /// Original thermal/material interface payload at the accepted endpoint.
    pub source: S,
    /// Existing numerical report, not a time-accuracy certificate.
    pub nonlinear: NewtonReport,
    /// Actual largest independently admitted cell residual [J].
    pub maximum_cell_residual_j: f64,
    /// All electrical evaluations, including discarded stencils and final check.
    pub circuit_trials: usize,
    /// Joules handed to the receiving cells, reconstructed from actual powers.
    pub deposited_joule_j: f64,
    /// Authored and callback source/boundary heat, EXCLUDING internal Joule heat.
    pub non_electrical_heat_j: f64,
    /// Driver storage + thermal storage + declared numerical dissipation
    /// - driver work - other external heat [J]. A field-motion driver includes
    /// body kinetic energy and mechanical work. Numerical dissipation is
    /// reported by its original driver step and is NEVER deposited as heat.
    pub total_energy_residual_j: f64,
}
/// Both endpoints belong to ONE physical interval; neither input is mutated.
pub struct CircuitEnthalpyCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    report: CircuitEnthalpyReport<S>,
}
impl<'m, S> CircuitEnthalpyCandidate<'m, S> {
    /// Proposed complete thermal history, including unchanged event cursor.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Proposed complete electrical history.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Source records and independently admitted local/global work.
    pub const fn report(&self) -> &CircuitEnthalpyReport<S> { &self.report }
    /// Publish both histories together; all electrical heat is already included.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, CircuitEnthalpyReport<S>) {
        (self.thermal, self.circuit, self.report)
    }
}

/// Solve one coupled multi-body electrical/enthalpy interval.
///
/// All resistors require a heat destination; passive thermal cells are allowed.
/// Resolve conductivity, resistivity and any extra sources against the SAME
/// proposed endpoint. Unknowns/initial guesses are cell-energy increments [J]
/// in sorted CellId order, not circuit branch order. Extra source maps must not
/// include Joule heat. Callbacks must not mutate external physical history.
///
/// Every circuit proposal starts from identical accepted charge/flux. Its
/// midpoint branch current produces heat using endpoint material resistance.
/// The existing thermal owner applies that heat, conduction and reservoirs in
/// its backward-Euler residual. Fresh local/global thermal and electrical work
/// checks run after Newton; no small-update or coefficient-lag fallback exists.
///
/// Geometry, circuit topology, C/L and caloric laws stay fixed. Process start
/// events first and split at interior events. Target events remain pending in
/// the accepted left-limit state. Source discontinuities require interval cuts.
/// Finite-difference stencils must remain supported by all material laws.
/// Numerical convergence does not prove uniqueness or physical stability.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: CircuitEnthalpyInput<'_>, controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<CircuitEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    driven_candidate(thermal, target_s, input, controls, initial_energy_increments_j, |state, _remaining| {
        let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
        let circuit = input.circuit.step_candidate(coefficients.resistors, input.ports,
            target_s-thermal.time_s(), controls.circuit, || (*cancel.borrow_mut())())?;
        Ok(DrivenCoefficients { circuit, sources_w:coefficients.sources_w,
            conductances_w_k:coefficients.conductances_w_k, source:coefficients.source, circuit_trials:1 })
    }, || (*cancel.borrow_mut())())
}

// A coupled driver may determine circuit ports from another physical owner.
// Keep branch heat deposition, cell residuals and the Newton/acceptance path
// identical to the ordinary circuit path. This is crate-private: callers cannot
// substitute a forged electrical endpoint or turn a trial into live history.
pub(crate) struct DrivenCoefficients<S, D = MaterialCircuitCandidate> {
    pub circuit: D,
    pub sources_w: BTreeMap<CellId, f64>,
    pub conductances_w_k: BTreeMap<u64, f64>,
    pub source: S,
    pub circuit_trials: usize,
}

// Both drivers expose their ORIGINAL accepted candidate and loss records. This
// private boundary only shares deposition/Newton/accounting; it cannot author
// another circuit integrator, change material laws or publish partial history.
pub(crate) struct DriverAccount<'a> {
    pub resistors: &'a [crate::material_circuit::ResistorHeat],
    pub contacts: &'a [crate::material_circuit::contacts::ContactHeat],
    pub heat_j: f64,
    /// Original time-discretization loss, separate from constitutive heating.
    pub numerical_dissipation_j: f64,
    pub stored_change_j: f64,
    pub supplied_j: f64,
}
pub(crate) trait ThermalDriver {
    type Model;
    type Step;
    fn account(&self) -> DriverAccount<'_>;
    fn into_parts(self) -> (Self::Model, Self::Step);
}
impl ThermalDriver for MaterialCircuitCandidate {
    type Model = MaterialCircuit;
    type Step = MaterialCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount { resistors: &step.resistors, contacts: &step.contacts,
            heat_j: step.heat_j, numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_after_j-step.stored_before_j,
            supplied_j: step.native.supplied }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        MaterialCircuitCandidate::into_parts(self)
    }
}
#[allow(clippy::too_many_arguments)]
pub(crate) fn driven_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: CircuitEnthalpyInput<'_>, controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>, usize) -> Result<DrivenCoefficients<S>, CircuitEnthalpyError<E>>,
    cancelled: impl FnMut() -> bool,
) -> Result<CircuitEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    driven_with_contacts_candidate(thermal, target_s, input, &[], controls,
        initial_energy_increments_j, resolve, cancelled)
}

// One physical/numerical owner for bulk-only and contact-heated networks.
#[allow(clippy::too_many_arguments)]
pub(crate) fn driven_with_contacts_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: CircuitEnthalpyInput<'_>, contact_cells: &[contacts::ContactHeatCells],
    controls: CircuitEnthalpySolve, initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>, usize) -> Result<DrivenCoefficients<S>, CircuitEnthalpyError<E>>,
    cancelled: impl FnMut() -> bool,
) -> Result<CircuitEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    let (thermal, circuit, report) = coupled_driver_candidate(thermal, target_s,
        input, contact_cells, controls, initial_energy_increments_j, resolve, cancelled)?;
    Ok(CircuitEnthalpyCandidate { thermal, circuit, report })
}

// The same residual and fresh acceptance for a fixed circuit or a circuit with
// its own reciprocally coupled body. Existing public circuit paths retain the
// exact loss deposition, controls, source callbacks and arithmetic below.
#[allow(clippy::too_many_arguments)]
pub(crate) fn coupled_driver_candidate<'m, S, E, D: ThermalDriver>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: CircuitEnthalpyInput<'_>, contact_cells: &[contacts::ContactHeatCells],
    controls: CircuitEnthalpySolve, initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>, usize) -> Result<DrivenCoefficients<S, D>, CircuitEnthalpyError<E>>,
    cancelled: impl FnMut() -> bool,
) -> Result<(DepositionSimulation<'m>, D::Model, CircuitEnthalpyReport<S, D::Step>), CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let checkpoint = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    checkpoint()?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let mut destinations = BTreeMap::new();
    if input.destinations.len() != input.circuit.resistors().len() {
        return Err(CircuitEnthalpyError::Invalid("every resistor requires exactly one thermal destination"));
    }
    for binding in input.destinations {
        checkpoint()?;
        if !thermal.cells().contains_key(&binding.cell) || destinations.insert(binding.branch, binding.cell).is_some()
            || !input.circuit.resistors().iter().any(|r| r.branch == binding.branch) {
            return Err(CircuitEnthalpyError::Invalid("duplicate, foreign or inactive resistor destination"));
        }
    }
    let contact_routes = contacts::admit_routes(thermal, input.circuit, contact_cells, checkpoint)?;
    cell_driver_candidate(thermal, target_s, controls, initial_energy_increments_j, resolve,
        |driver: &D, dt, check| {
        let account = driver.account();
        let mut powers = BTreeMap::new();
        for branch in account.resistors {
            check()?;
            let cell = *destinations.get(&branch.material.branch)
                .ok_or(CircuitEnthalpyError::Invalid("circuit returned heat without a destination"))?;
            let power = finite(branch.heat_j/dt)?;
            if branch.heat_j != 0.0 && power == 0.0 {
                return Err(CircuitEnthalpyError::Invalid("branch heat disappeared during power conversion"));
            }
            let sum = powers.entry(cell).or_insert(0.0);
            *sum = finite(*sum + power)?;
        }
        contacts::deposit(account.contacts, &contact_routes, dt, &mut powers, || check())?;
        Ok(powers)
    }, || (*cancel.borrow_mut())())
}

fn admit_driver_unknowns<E>(
    thermal: &DepositionSimulation<'_>, controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<(), CircuitEnthalpyError<E>> {
    controls.newton.validate()?;
    if controls.maximum_cells == 0 || thermal.cells().len() > controls.maximum_cells
        || initial_energy_increments_j.len() != thermal.cells().len()
        || controls.maximum_circuit_trials == 0 || controls.maximum_newton_iterations == 0
        || [controls.cell_energy_tolerance_j, controls.thermal_energy_tolerance_j, controls.total_energy_tolerance_j]
            .iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(CircuitEnthalpyError::Invalid("invalid thermal unknowns, work budget or absolute energy tolerance"));
    }
    Ok(())
}

// Only physical heat routing differs between branch/contact and spatial winding
// drivers. The thermal residual, numerical owner, budget and fresh acceptance
// below remain shared; this boundary is private and accepts admitted candidates.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cell_driver_candidate<'m, S, E, D: ThermalDriver>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        usize,
    ) -> Result<DrivenCoefficients<S, D>, CircuitEnthalpyError<E>>,
    deposit: impl Fn(
        &D,
        f64,
        &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<E>>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (
        DepositionSimulation<'m>,
        D::Model,
        CircuitEnthalpyReport<S, D::Step>,
    ),
    CircuitEnthalpyError<E>,
> {
    let cancel = RefCell::new(cancelled);
    let checkpoint = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    checkpoint()?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
    let dt = step.duration_s();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| step.evaluate(x, |state| {
        checkpoint()?;
        if trials.get() >= controls.maximum_circuit_trials { return Err(CircuitEnthalpyError::CircuitBudget); }
        let remaining = controls.maximum_circuit_trials-trials.get();
        let coefficients = resolve(state, remaining)?;
        if coefficients.circuit_trials == 0 || coefficients.circuit_trials > remaining {
            return Err(CircuitEnthalpyError::CircuitBudget);
        }
        trials.set(trials.get()+coefficients.circuit_trials);
        let circuit = coefficients.circuit;
        let account = circuit.account();
        if !account.numerical_dissipation_j.is_finite() || account.numerical_dissipation_j < 0.0 {
            return Err(CircuitEnthalpyError::Invalid("driver numerical dissipation must be finite and nonnegative"));
        }
        let mut powers = deposit(&circuit, dt, &mut || checkpoint())?;
        let deposited = finite(dt * powers.values().try_fold(0.0, |s, &p| finite(s+p))?)?;
        if (deposited-account.heat_j).abs() > controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid("branch heat and receiving-cell power disagree"));
        }
        if coefficients.sources_w.len() > state.cells().len() {
            return Err(CircuitEnthalpyError::Invalid("too many external source cells"));
        }
        for (cell, power) in coefficients.sources_w {
            checkpoint()?; finite(power)?;
            if !state.cells().contains_key(&cell) { return Err(CircuitEnthalpyError::Invalid("inactive external source cell")); }
            let sum = powers.entry(cell).or_insert(0.0);
            *sum = finite(*sum + power)?;
        }
        Ok::<_, CircuitEnthalpyError<E>>(StepCoefficients { sources_w: powers,
            conductances_w_k: coefficients.conductances_w_k,
            payload: (circuit, coefficients.source, deposited) })
    }, checkpoint);
    let problem = FallibleFiniteDifference::new(step.cell_ids().len(), controls.difference_step_j,
        |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec())).map_err(numerical)?;
    let construction = NewtonKrylovState::new(&problem, initial_energy_increments_j.to_vec(), controls.newton);
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
        &mut || (*cancel.borrow_mut())() || problem.has_failed());
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete { return Err(CircuitEnthalpyError::Interrupted(result.progress)); }
    if !result.report.converged { return Err(CircuitEnthalpyError::NotConverged(result.report)); }
    let trial = evaluate(&solver.x)?;
    let maximum_cell_residual_j = trial.maximum_cell_residual_j();
    let (next, (candidate, source, deposited_joule_j), audit) = trial.accept(
        controls.cell_energy_tolerance_j, controls.thermal_energy_tolerance_j, checkpoint)?;
    let non_electrical_heat_j = finite(audit.external_heat_j - deposited_joule_j)?;
    let account = candidate.account();
    let total = finite(account.stored_change_j + account.numerical_dissipation_j
        + audit.stored_energy_change_j-account.supplied_j-non_electrical_heat_j)?;
    if total.abs() > controls.total_energy_tolerance_j {
        return Err(CircuitEnthalpyError::Invalid("circuit and thermal storage do not close the supplied energy"));
    }
    checkpoint()?;
    let (circuit, electrical) = candidate.into_parts();
    Ok((next, circuit, CircuitEnthalpyReport { electrical, thermal: audit, source, nonlinear: result.report,
            maximum_cell_residual_j, circuit_trials: trials.get(), deposited_joule_j,
            non_electrical_heat_j, total_energy_residual_j: total }))
}
