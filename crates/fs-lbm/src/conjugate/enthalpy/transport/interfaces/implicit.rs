//! Closed transient electrothermal feedback, not a one-pass coefficient lag.
//!
//! The existing Newton--Krylov owner solves the shared thermal owner's
//! backward-Euler residual in cell energy increments [J]. Every evaluation
//! re-solves the existing electrical/contact problem at the SAME trial h/T,
//! including phase change at constant Kelvin. No trial advances accepted
//! history or deposits its Joule work. Only the converged endpoint's power is
//! integrated, once, over this interval.
//!
//! Jacobian actions are explicitly numerical directional differences in joule
//! coordinates, not analytic material tangents or adjoints. The caller selects
//! the difference scale, solver controls, initial guess and work budgets.
//! Unsupported stencils refuse; neither source extrapolation nor an explicit
//! fallback is used. Smoothness, uniqueness, nonlinear convergence and temporal
//! accuracy are not inferred from algebraic convergence.

use super::*;
use std::cell::{Cell, RefCell};

/// Error-controlled evolution using the existing implicit step and energy accounts.
pub mod adaptive;
use fs_conduction::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};

/// Separate nonlinear work, finite-difference scale and physical residual gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitJouleConfig {
    /// Maximum active thermal cells, including cells outside the current map.
    pub max_unknowns: usize,
    /// Bounded outer Newton attempts; inner work is bounded by `newton`.
    pub max_newton_iterations: usize,
    /// Maximum-component perturbation for a numerical Jacobian action [J].
    pub difference_step_j: f64,
    /// Absolute maximum per-cell backward-Euler energy defect [J].
    pub cell_energy_tolerance_j: f64,
    /// Absolute total thermal-energy defect [J], independent of electrical error.
    pub total_energy_tolerance_j: f64,
    /// Existing solver policy. Its residual coordinates are joules.
    pub newton: NewtonKrylovConfig,
}

/// Original physical refusals are retained even inside infallible solver hooks.
#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitJouleError {
    /// Original material, electrical, thermal or current-budget error.
    Physical(JouleEnthalpyError),
    /// Invalid numerical coupling controls.
    Invalid(&'static str),
    /// Existing Newton construction/configuration error.
    Newton(NewtonError),
    /// Existing solver did not converge; no thermal endpoint is returned.
    NotConverged(NewtonReport),
    /// Interrupted solver work or a contained numerical callback panic.
    Interrupted(SolverRunProgress),
}
impl From<JouleEnthalpyError> for ImplicitJouleError {
    fn from(e: JouleEnthalpyError) -> Self { Self::Physical(e) }
}
impl From<DepositionError> for ImplicitJouleError {
    fn from(e: DepositionError) -> Self { Self::Physical(e.into()) }
}
impl From<NewtonError> for ImplicitJouleError {
    fn from(e: NewtonError) -> Self { Self::Newton(e) }
}
impl core::fmt::Display for ImplicitJouleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Physical(e) => write!(f, "implicit electrothermal feedback: {e}"),
            Self::Invalid(s) => write!(f, "implicit electrothermal controls: {s}"),
            Self::Newton(e) => write!(f, "implicit electrothermal Newton: {e}"),
            Self::NotConverged(r) => write!(f, "implicit electrothermal residual {} J after {} attempts: {:?}",
                r.residual_norm, r.iterations, r.diagnosis),
            Self::Interrupted(p) => write!(f, "implicit electrothermal solve interrupted: {p:?}"),
        }
    }
}
impl std::error::Error for ImplicitJouleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Physical(e) => Some(e), Self::Newton(e) => Some(e), _ => None }
    }
}

/// One accepted implicit interval; the endpoint fields actually supplied heat.
#[derive(Debug)]
pub struct ImplicitJouleReport<F, I, J> {
    /// Total heat includes authored, additional and Joule input exactly once.
    pub thermal: EnergyAudit,
    /// Existing Newton report, not a physical or timestep-accuracy certificate.
    pub nonlinear: NewtonReport,
    /// Fresh maximum local residual that passed the independent acceptance gate.
    pub maximum_cell_residual_j: f64,
    /// All current solves, including numerical Jacobian probes and final check.
    pub current_solves: usize,
    /// Actual implicit endpoint electrical state, not a start-state diagnostic.
    pub endpoint: EnthalpyCurrentState,
    /// Exact endpoint material/phase evaluations used for current and transport.
    pub material: Arc<EnthalpyTransportEpoch<F, I>>,
    /// Exact ordered interfaces used by that same current solve.
    pub interfaces: Arc<ElectricalInterfaceEpoch<J>>,
    /// Backward-Euler electrode work at the accepted endpoint [J].
    pub electrode_energy_j: f64,
    /// Irreversible electrical heat integrated over the same interval [J].
    pub joule_energy_j: f64,
    /// Joule heat actually handed to receiving cells [J].
    pub deposited_joule_j: f64,
}

// Share numerical stencils and source-error latching with circuit coupling.
fn numerical_error(error: FiniteDifferenceError<JouleEnthalpyError>) -> ImplicitJouleError {
    match error {
        FiniteDifferenceError::Source(e) => e.into(),
        FiniteDifferenceError::Invalid(reason) => ImplicitJouleError::Invalid(reason),
    }
}

/// Solve ONE closed implicit interval with actual bulk/interface/heat feedback.
///
/// `initial_energy_increments_j` is an initial guess in the thermal step's stable
/// CellId order, not voxel order. It is not an applied source. `resolve` and both
/// contact callbacks have the existing material-interface path's semantics.
/// `external_sources` adds independent endpoint powers, evaluated against that
/// same state, without advancing another thermal history. All callbacks must be
/// pure with respect to physical history; mutable observer state is not admitted.
///
/// Current-solve limits in `input` bound ALL nonlinear material epochs. A final
/// fresh endpoint evaluation and local/global energy gates are mandatory even
/// if the Newton stopping rule declares success. No returned candidate on late
/// failure, exhausted work, unsupported stencil, cancellation or nonconvergence.
///
/// Authored geometry, reference mass, boundaries, caloric laws and topology stay
/// fixed. Events inside the interval refuse; target-time events remain pending
/// until explicitly processed after coupled acceptance. Electrical drive/source
/// discontinuities require interval boundaries. This is stationary Ohmic and
/// control-volume backward Euler, not melt motion, induction or strong coupling
/// to unrepresented mechanical/history variables. Numerical Jacobian and time
/// refinement must be checked independently for the supplied material laws.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_joule_enthalpy_implicit<'t, F, I, J>(
    thermal: &DepositionSimulation<'t>, target_s: f64,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate,
    config: ImplicitJouleConfig, initial_energy_increments_j: &[f64],
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64,I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<ElectricalInterfaceResponse<J>, JouleEnthalpyError>,
    external_sources: impl Fn(&DepositionSimulation<'t>) -> Result<BTreeMap<CellId,f64>, JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>, ImplicitJouleReport<F,I,J>), ImplicitJouleError> {
    step_counted(thermal, target_s, input, link_ids, gate, config, initial_energy_increments_j,
        resolve, resolve_link, resolve_interface, external_sources, &Cell::new(0))
}

// Retain actually spent current work on every return, including failed Newton
// trials. Adaptive retries share one budget; no failed source becomes free work.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn step_counted<'t, F, I, J>(
    thermal: &DepositionSimulation<'t>, target_s: f64,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate,
    config: ImplicitJouleConfig, initial_energy_increments_j: &[f64],
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64,I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<ElectricalInterfaceResponse<J>, JouleEnthalpyError>,
    external_sources: impl Fn(&DepositionSimulation<'t>) -> Result<BTreeMap<CellId,f64>, JouleEnthalpyError>,
    current_work: &Cell<usize>,
) -> Result<(DepositionSimulation<'t>, ImplicitJouleReport<F,I,J>), ImplicitJouleError> {
    checkpoint(gate)?;
    config.newton.validate()?;
    if config.max_unknowns == 0 || thermal.cells().len() > config.max_unknowns
        || initial_energy_increments_j.len() != thermal.cells().len()
        || config.max_newton_iterations == 0 || !config.difference_step_j.is_finite() || config.difference_step_j <= 0.0
        || !config.cell_energy_tolerance_j.is_finite() || config.cell_energy_tolerance_j < 0.0
        || !config.total_energy_tolerance_j.is_finite() || config.total_energy_tolerance_j < 0.0 {
        return Err(ImplicitJouleError::Invalid("invalid unknown/work budget, difference scale or joule tolerance"));
    }
    admit_map(thermal, input, gate)?;
    admit_interfaces(input, gate)?;
    let step = BackwardEulerStep::new(thermal, target_s, ||checkpoint(gate))?;
    let source = RefCell::new(InterfaceSource {
        transport: TransportSource { account: Account::new(input,gate), selected:link_ids,
            indices: input.cell_ids.iter().copied().enumerate().map(|(c,id)|(id,c)).collect(),
            bulk: &resolve, link: &resolve_link, last: None },
        interface: &resolve_interface, last: None,
    });
    let evaluate = |x: &[f64]| step.evaluate(x, |state| {
        let evaluated = source.borrow_mut().evaluate(state);
        current_work.set(source.borrow().transport.account.solves);
        let (current, bulk, interfaces) = evaluated?;
        let mut sources_w: BTreeMap<_,_> = input.cell_ids.iter().copied().zip(current.electric.joule_w.iter().copied()).collect();
        let external = external_sources(state)?;
        if external.len() > state.cells().len() { return Err(invalid("too many implicit external sources")); }
        for (i,(id,q)) in external.into_iter().enumerate() {
            if i % 256 == 0 { checkpoint(gate)?; }
            finite(q)?;
            if !state.cells().contains_key(&id) { return Err(invalid("inactive implicit external-source cell")); }
            let total = sources_w.entry(id).or_insert(0.0);
            *total = finite(*total + q)?;
        }
        let conductances_w_k = bulk.links.iter().map(|s|(s.patch.link.id,s.conductance_w_k)).collect();
        Ok::<_,JouleEnthalpyError>(StepCoefficients { sources_w, conductances_w_k, payload: (current,bulk,interfaces) })
    }, ||checkpoint(gate));
    let problem = FallibleFiniteDifference::new(step.cell_ids().len(), config.difference_step_j,
        |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec())).map_err(numerical_error)?;
    let construction = NewtonKrylovState::new(&problem, initial_energy_increments_j.to_vec(), config.newton);
    if let Some(error) = problem.take_failure() { return Err(numerical_error(error)); }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, config.max_newton_iterations,
        &mut || gate.is_requested() || problem.has_failed());
    if let Some(error) = problem.take_failure() { return Err(numerical_error(error)); }
    checkpoint(gate)?;
    if result.progress != SolverRunProgress::Complete { return Err(ImplicitJouleError::Interrupted(result.progress)); }
    if !result.report.converged { return Err(ImplicitJouleError::NotConverged(result.report)); }
    let trial = evaluate(&solver.x)?;
    let maximum_cell_residual_j = trial.maximum_cell_residual_j();
    let (next, (endpoint, material, interfaces), audit) = trial.accept(
        config.cell_energy_tolerance_j, config.total_energy_tolerance_j, ||checkpoint(gate))?;
    let dt = step.duration_s();
    let electrode_w = endpoint.electric.electrodes.iter().try_fold(0.0,
        |s,e| finite(s + electrode_power(e)?))?;
    let deposited_w = endpoint.electric.joule_w.iter().try_fold(0.0, |s,&q|finite(s+q))?;
    let electrode_energy_j = finite(dt*electrode_w)?;
    let joule_energy_j = finite(dt*endpoint.electric.total_joule_w)?;
    let deposited_joule_j = finite(dt*deposited_w)?;
    if (electrode_energy_j-joule_energy_j).abs() > input.energy_tolerance_j
        || (joule_energy_j-deposited_joule_j).abs() > input.energy_tolerance_j {
        return Err(invalid("implicit electrical work and deposited heat do not close").into());
    }
    checkpoint(gate)?;
    let current_solves = source.borrow().transport.account.solves;
    Ok((next, ImplicitJouleReport { thermal:audit, nonlinear:result.report,
        maximum_cell_residual_j, current_solves, endpoint, material, interfaces,
        electrode_energy_j, joule_energy_j, deposited_joule_j }))
}
