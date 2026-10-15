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
use std::cell::RefCell;
use fs_conduction::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    NonlinearProblem, SolverRunProgress};

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

// The nonlinear owner has infallible numerical callbacks. Retain the original
// typed refusal out of band and poison all subsequent numerical work, rather
// than returning a fabricated zero residual or replacing unsupported material.
struct NumericalResidual<E> {
    evaluate: E,
    n: usize,
    difference_j: f64,
    failure: RefCell<Option<JouleEnthalpyError>>,
}
impl<E: Fn(&[f64]) -> Result<Vec<f64>, JouleEnthalpyError>> NumericalResidual<E> {
    fn checked(&self, x: &[f64]) -> Result<Vec<f64>, JouleEnthalpyError> {
        if let Some(e) = self.failure.borrow().as_ref() { return Err(e.clone()); }
        (self.evaluate)(x)
    }
    fn refuse(&self, error: JouleEnthalpyError, output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.failure.borrow().is_none() { *self.failure.borrow_mut() = Some(error); }
    }
    fn take_failure(&self) -> Result<(), ImplicitJouleError> {
        match self.failure.borrow_mut().take() { Some(e) => Err(e.into()), None => Ok(()) }
    }
}
impl<E: Fn(&[f64]) -> Result<Vec<f64>, JouleEnthalpyError>> NonlinearProblem for NumericalResidual<E> {
    fn dimension(&self) -> usize { self.n }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        match self.checked(x) {
            Ok(values) => output.copy_from_slice(&values),
            Err(e) => self.refuse(e, output),
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        let result = (|| {
            let scale = direction.iter().copied().map(f64::abs).fold(0.0, f64::max);
            if !direction.iter().all(|x| x.is_finite()) { return Err(invalid("nonfinite Jacobian direction")); }
            if scale == 0.0 { return Ok(vec![0.0; self.n]); }
            let delta = self.difference_j / scale;
            if !delta.is_finite() || delta <= 0.0 { return Err(invalid("unrepresentable Jacobian difference scale")); }
            let shifted: Vec<_> = x.iter().zip(direction).map(|(&a,&d)| a + delta*d).collect();
            if !shifted.iter().all(|a| a.is_finite()) || shifted.as_slice() == x {
                return Err(invalid("Jacobian stencil is nonfinite or too small to change energy"));
            }
            let base = self.checked(x)?;
            let other = self.checked(&shifted)?;
            base.iter().zip(&other).map(|(&a,&b)| finite((b-a)/delta)).collect()
        })();
        match result {
            Ok(values) => output.copy_from_slice(&values),
            Err(e) => self.refuse(e, output),
        }
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
        let (current, bulk, interfaces) = source.borrow_mut().evaluate(state)?;
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
    let problem = NumericalResidual { evaluate: |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec()),
        n: step.cell_ids().len(), difference_j: config.difference_step_j, failure: RefCell::new(None) };
    let construction = NewtonKrylovState::new(&problem, initial_energy_increments_j.to_vec(), config.newton);
    problem.take_failure()?;
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, config.max_newton_iterations,
        &mut || gate.is_requested() || problem.failure.borrow().is_some());
    problem.take_failure()?;
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
