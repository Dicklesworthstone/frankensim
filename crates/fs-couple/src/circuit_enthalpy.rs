//! Transient Kirchhoff circuits coupled to the shared enthalpy network.
//!
//! Every material resistor has an explicit receiving thermal cell. Conduction,
//! reservoirs, latent storage and local energy residuals stay with fs-thermal;
//! charges, fluxes, KCL/KVL and branch losses stay with MaterialCircuit. No trial
//! advances either accepted history. This is backward Euler for thermal rates
//! with circuit midpoint currents and ENDPOINT material coefficients: first-order
//! coupled accuracy in general, not a second-order or unconditional-stability claim.

/// Data-only caloric/card, resistor and conductive-link bindings.
pub mod material;

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
pub struct CircuitEnthalpyReport<S> {
    /// Original circuit step: finite storage, branch heat and external port work.
    pub electrical: MaterialCircuitStep,
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
    /// Circuit storage + thermal storage - port work - non-electrical heat [J].
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
    let checkpoint = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    checkpoint()?;
    controls.newton.validate()?;
    if controls.maximum_cells == 0 || thermal.cells().len() > controls.maximum_cells
        || initial_energy_increments_j.len() != thermal.cells().len()
        || controls.maximum_circuit_trials == 0 || controls.maximum_newton_iterations == 0
        || [controls.cell_energy_tolerance_j, controls.thermal_energy_tolerance_j, controls.total_energy_tolerance_j]
            .iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(CircuitEnthalpyError::Invalid("invalid thermal unknowns, work budget or absolute energy tolerance"));
    }
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
    let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
    let dt = step.duration_s();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| step.evaluate(x, |state| {
        checkpoint()?;
        if trials.get() >= controls.maximum_circuit_trials { return Err(CircuitEnthalpyError::CircuitBudget); }
        trials.set(trials.get()+1);
        let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
        let circuit: MaterialCircuitCandidate = input.circuit.step_candidate(coefficients.resistors,
            input.ports, dt, controls.circuit, || (*cancel.borrow_mut())())?;
        let mut powers = BTreeMap::new();
        for branch in &circuit.step().resistors {
            checkpoint()?;
            let cell = *destinations.get(&branch.material.branch)
                .ok_or(CircuitEnthalpyError::Invalid("circuit returned heat without a destination"))?;
            let power = finite(branch.heat_j/dt)?;
            if branch.heat_j != 0.0 && power == 0.0 {
                return Err(CircuitEnthalpyError::Invalid("branch heat disappeared during power conversion"));
            }
            let sum = powers.entry(cell).or_insert(0.0);
            *sum = finite(*sum + power)?;
        }
        let deposited = finite(dt * powers.values().try_fold(0.0, |s, &p| finite(s+p))?)?;
        if (deposited-circuit.step().heat_j).abs() > controls.total_energy_tolerance_j {
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
    let (circuit, electrical) = candidate.into_parts();
    let non_electrical_heat_j = finite(audit.external_heat_j - deposited_joule_j)?;
    let total = finite((electrical.stored_after_j-electrical.stored_before_j)
        + audit.stored_energy_change_j-electrical.native.supplied-non_electrical_heat_j)?;
    if total.abs() > controls.total_energy_tolerance_j {
        return Err(CircuitEnthalpyError::Invalid("circuit and thermal storage do not close the supplied energy"));
    }
    checkpoint()?;
    Ok(CircuitEnthalpyCandidate { thermal: next, circuit,
        report: CircuitEnthalpyReport { electrical, thermal: audit, source, nonlinear: result.report,
            maximum_cell_residual_j, circuit_trials: trials.get(), deposited_joule_j,
            non_electrical_heat_j, total_energy_residual_j: total } })
}
