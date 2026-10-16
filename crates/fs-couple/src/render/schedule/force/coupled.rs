//! Two-way elastic/viscous connections between retained modal components.
//!
//! Every retained mode, bilateral spring and dashpot enters one complete
//! quadratic generator, advanced by exact held-external-force propagation.
//! Independent integrated resistance channels retain component and connection
//! loss accounting. Nonlinear contact/friction wraps this same affine response
//! with its existing held effort times displacement work. Static preload keeps
//! the existing small Schur solve. Dense preparation has explicit work/capacity
//! refusals; no modes are truncated. Pressure transfers retain their original
//! narrow-band, read-only scope; radiation loading is not added.

use fs_exec::CancelGate;
use fs_la::factor::{Cholesky, cholesky};

use crate::modal_acoustic_time::{
    ModalAcousticState, ModalAcousticTimeError, ModalAcousticTimeModel,
};

mod exact;

/// One attachment in a component's unchanged mass-normalized modal basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalAttachment {
    /// Component index, in construction order.
    pub component: usize,
    /// Displacement/force participation [1/sqrt(kg)], one entry per mode.
    pub shapes: Vec<f64>,
}

/// A bilateral spring and dashpot between two explicitly declared attachments.
/// Its extension is left displacement minus right displacement minus rest.
/// This is not unilateral contact, friction, a rigid constraint or a preset.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalConnection {
    /// Positive side of the relative-displacement coordinate.
    pub left: ModalAttachment,
    /// Negative side. Zero shapes explicitly describe a fixed support.
    pub right: ModalAttachment,
    /// Nonnegative spring stiffness [N/m]. Zero leaves only the dashpot.
    pub stiffness_n_m: f64,
    /// Nonnegative viscous resistance [N s/m]. Zero is an elastic connection.
    pub damping_n_s_m: f64,
    /// Stress-free relative displacement [m]. May be signed.
    pub rest_extension_m: f64,
}

/// One bilateral spring/dashpot on a complete signed modal coordinate.
/// Its extension is `column^T q - rest_extension_m`; the SAME column maps
/// its scalar reaction back to every participating body. Thus a shared cavity
/// or distributed damper retains all cross terms without merging components
/// or replacing one physical port with pairwise connections.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalColumnConnection {
    /// Displacement/force participation [1/sqrt(kg)], in original component
    /// order followed by each component's mode order. Include exact zeros for
    /// unattached modes. The column is not normalized or regrouped.
    pub column: Vec<f64>,
    /// Nonnegative spring stiffness [N/m]. Zero leaves only the dashpot.
    pub stiffness_n_m: f64,
    /// Nonnegative viscous resistance [N s/m]. Zero is an elastic connection.
    pub damping_n_s_m: f64,
    /// Stress-free displacement of the complete signed coordinate [m].
    pub rest_extension_m: f64,
}

struct ConnectionLaw {
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    rest_extension_m: f64,
}

/// Explicit numerical and work admission for a fixed coupled configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalCouplingConfig {
    /// Maximum total component modes (hard ceiling 4096).
    pub max_modes: usize,
    /// Maximum connections (hard ceiling 64).
    pub max_connections: usize,
    /// Budget for n_modes * (n_connections+1)^2 setup terms.
    pub max_setup_terms: usize,
    /// Guard applied to an algebraic upper bound on coupled natural frequency.
    /// Ordinary floating-point evaluation, not an interval certificate.
    pub nyquist_guard_fraction: f64,
    /// Absolute combined modal PLUS connection-spring energy ceiling [J].
    pub maximum_total_energy_j: f64,
    /// Absolute pressure ceiling on the sum of component observations [Pa].
    pub maximum_abs_pressure_pa: f64,
    /// Absolute reaction-force ceiling on each connection [N].
    pub maximum_abs_connection_force_n: f64,
    /// Scaled infinity-norm residual gate for the connection solve, in (0,1).
    pub solve_relative_tolerance: f64,
    /// Absolute discrete energy-closure tolerance [J].
    pub energy_absolute_tolerance_j: f64,
    /// Relative discrete energy-closure tolerance, in [0,1).
    pub energy_relative_tolerance: f64,
}

/// A refused trial never changes accepted component states or the sample clock.
#[derive(Debug)]
pub enum ModalCouplingError {
    /// Complete coupled linear preparation or propagation refused.
    Propagation(fs_phs::PhsError),
    /// Invalid shapes, parameters, derived numbers, or work admission.
    Invalid { what: &'static str },
    /// Cancellation observed before publishing a complete coupled step.
    Cancelled,
    /// An original component rejected its candidate state or transition.
    Component { component: usize, source: ModalAcousticTimeError },
    /// fs-la could not factor the finite symmetric connection matrix.
    Factor(fs_la::factor::FactorError),
    /// The existing unilateral-contact law refused a trial.
    ContactLaw(fs_dcontact::DContactError),
    /// A contact root or its final constitutive check exhausted the declared tolerance.
    ContactSolve { residual_n: f64, tolerance_n: f64, iterations: usize },
    /// Recomputed connection residual exceeds the caller's tolerance.
    SolveResidual { relative: f64, tolerance: f64 },
    /// A combined physical ceiling was exceeded.
    Budget { what: &'static str, value: f64, limit: f64 },
    /// Storage change plus all disclosed losses does not match external work.
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}

impl core::fmt::Display for ModalCouplingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Propagation(error) => write!(f,"coupled linear propagation: {error}"),
            Self::Invalid { what } => write!(f, "modal coupling: {what}"),
            Self::Cancelled => write!(f, "modal coupling cancelled before publication"),
            Self::Component { component, source } => write!(f, "modal component {component}: {source}"),
            Self::ContactLaw(error) => write!(f, "modal contact law: {error}"),
            Self::ContactSolve { residual_n, tolerance_n, iterations } => write!(f,
                "modal contact residual {residual_n:e} N exceeds {tolerance_n:e} N after {iterations} iterations"),
            Self::Factor(error) => write!(f, "modal connection factor: {error}"),
            Self::SolveResidual { relative, tolerance } => write!(f, "connection residual {relative:e} exceeds {tolerance:e}"),
            Self::Budget { what, value, limit } => write!(f, "coupled {what}: {value:e} exceeds {limit:e}"),
            Self::EnergyBalance { residual_j, tolerance_j } => write!(f, "coupled energy residual {residual_j:e} J exceeds {tolerance_j:e} J"),
        }
    }
}
impl std::error::Error for ModalCouplingError {}

/// Diagnostics for the last accepted sample, never a partially staged trial.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledModalFrame {
    /// One-based accepted sample ordinal.
    pub sample: u64,
    /// Sum of the original component pressure observations [Pa].
    pub observer_pressure_pa: f64,
    /// Energy in retained oscillators, excluding connection springs [J].
    pub modal_energy_j: f64,
    /// Energy in connection springs [J].
    pub connection_energy_j: f64,
    /// Work of the externally held generalized forces only [J].
    pub external_work_j: f64,
    /// Independent time integral of all original component viscous losses [J].
    pub component_dissipation_j: f64,
    /// Independent time integral of c*(extension velocity)^2 [J].
    pub connection_dissipation_j: f64,
    /// Combined storage change + dissipation - external work [J].
    pub energy_residual_j: f64,
    /// Actual absolute + relative energy tolerance used by the gate [J].
    pub energy_tolerance_j: f64,
    /// Connection-equation residual for legacy accounts. Exact full-pencil
    /// propagation has no connection solve and records zero here; its actual
    /// numerical check is the independently integrated energy residual above.
    pub solve_relative_residual: f64,
    /// Endpoint signed reaction multiplying each connection's signed column [N].
    /// For a pair connection this is the force on its LEFT attachment.
    pub connection_forces_n: Vec<f64>,
}

/// A fixed network with two-way reactions and transactional accepted states.
/// Models, their individual budgets, damping and pressure transfers are retained.
/// Exact dense propagation allocates bounded matrix scratch during stepping.
/// This is not a measured real-time or experimental-validity claim.
pub struct CoupledModalSystem {
    models: Vec<ModalAcousticTimeModel>,
    candidates: Vec<ModalAcousticTimeModel>,
    exact: Option<exact::ExactNetwork>,
    offsets: Vec<usize>,
    connections: Vec<ConnectionLaw>,
    columns: Vec<Vec<f64>>,
    config: ModalCouplingConfig,
    dt: f64,
    old_q: Vec<f64>,
    free_delta: Vec<f64>,
    forces: Vec<f64>,
    reactions: Vec<f64>,
    frame: CoupledModalFrame,
    staged_frame: CoupledModalFrame,
    ports: crate::port_adoption::ModulePorts,
    staged_ports: crate::port_adoption::ModulePorts,
}

impl CoupledModalSystem {
    /// Prepare the small connection solve over caller-owned modal components.
    /// Restored component states are kept; no equilibrium or tuning is invented.
    /// A conservative algebraic stiffness bound rejects potentially aliased
    /// coupled modes. It may refuse a configuration a tighter eigensolve admits.
    pub fn new(
        models: Vec<ModalAcousticTimeModel>,
        connections: Vec<ModalConnection>,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let (offsets, count, _) = connection_layout(&models, connections.len(), config, gate)?;
        let mut complete = Vec::with_capacity(connections.len());
        for link in connections {
            poll(Some(gate))?;
            let mut column = vec![0.0; count];
            for (attachment, sign) in [(&link.left, 1.0), (&link.right, -1.0)] {
                let model = models.get(attachment.component).ok_or_else(|| invalid("attachment names an unknown component"))?;
                if attachment.shapes.len() != model.modes().len() || attachment.shapes.iter().any(|x| !x.is_finite()) {
                    return Err(invalid("attachment must match its component's finite mass-normalized shapes"));
                }
                for (mode, shape) in attachment.shapes.iter().enumerate() {
                    let index = offsets[attachment.component] + mode;
                    column[index] = finite(column[index] + sign * shape)?;
                }
            }
            complete.push(ModalColumnConnection {
                column, stiffness_n_m: link.stiffness_n_m,
                damping_n_s_m: link.damping_n_s_m, rest_extension_m: link.rest_extension_m,
            });
        }
        Self::new_with_columns(models, complete, config, gate)
    }

    /// Prepare simultaneous bilateral reactions on complete signed columns.
    /// Each column describes ONE physical displacement coordinate, even when
    /// several original bodies participate. Component states, damping, pressure
    /// transfers and individual budgets retain their original ownership.
    /// The same admission and finite-step solve as [`Self::new`] apply.
    pub fn new_with_columns(
        models: Vec<ModalAcousticTimeModel>,
        connections: Vec<ModalColumnConnection>,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let (offsets, count, dt) = connection_layout(&models, connections.len(), config, gate)?;
        let links = connections.len();
        let mut columns = Vec::with_capacity(links);
        let mut laws = Vec::with_capacity(links);
        for link in connections {
            poll(Some(gate))?;
            if !link.stiffness_n_m.is_finite() || link.stiffness_n_m < 0.0
                || !link.damping_n_s_m.is_finite() || link.damping_n_s_m < 0.0
                || !link.rest_extension_m.is_finite() {
                return Err(invalid("connection stiffness/damping must be finite nonnegative; rest must be finite"));
            }
            if link.column.len() != count || link.column.iter().any(|b| !b.is_finite()) {
                return Err(invalid("connection column must match every original mode with finite signed participation"));
            }
            columns.push(link.column);
            laws.push(ConnectionLaw { stiffness_n_m: link.stiffness_n_m,
                damping_n_s_m: link.damping_n_s_m, rest_extension_m: link.rest_extension_m });
        }
        let connections = laws;
        let mut upper_squared = 0.0_f64;
        for model in &models {
            poll(Some(gate))?;
            for &mode in model.modes() {
                upper_squared = upper_squared.max(finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s)?);
            }
        }
        // lambda_max(diag(omega^2) + sum k*b*b^T) <= max(omega^2) + sum k*||b||^2.
        // This is an admission screen, not a rounded interval eigenvalue proof.
        for (column, link) in columns.iter().zip(&connections) {
            let mut squared_norm = 0.0;
            for &b in column { squared_norm = finite(squared_norm + b*b)?; }
            upper_squared = finite(upper_squared + link.stiffness_n_m * squared_norm)?;
        }
        let guard = core::f64::consts::PI / dt * config.nyquist_guard_fraction;
        if upper_squared.sqrt() > guard { return Err(invalid("coupled stiffness bound exceeds the Nyquist guard")); }
        let exact = if connections
            .iter()
            .any(|link| link.stiffness_n_m != 0.0 || link.damping_n_s_m != 0.0)
        {
            Some(exact::ExactNetwork::new(
                &models,
                &columns,
                &connections,
                count,
                gate,
            )?)
        } else {
            None
        };
        let candidates = models.clone();
        let frame = CoupledModalFrame {
            sample: 0, observer_pressure_pa: 0.0, modal_energy_j: 0.0, connection_energy_j: 0.0,
            external_work_j: 0.0, component_dissipation_j: 0.0, connection_dissipation_j: 0.0,
            energy_residual_j: 0.0, energy_tolerance_j: 0.0, solve_relative_residual: 0.0,
            connection_forces_n: vec![0.0; links],
        };
        let mut ports = crate::port_adoption::ModulePorts::mechanical("coupled-modal", links);
        // Each signed B_i^T q_i is an explicit metre-valued work chart.
        // Add the connection-storage endpoint to the participating bodies.
        for (j, column) in columns.iter().enumerate() {
            let bodies = offsets
                .windows(2)
                .filter(|r| column[r[0]..r[1]].iter().any(|&b| b != 0.0))
                .count();
            if bodies != 0 {
                ports.junctions.push(crate::port_adoption::junction(
                    "coupled-modal",
                    &format!("column/{j}"),
                    bodies + 1,
                    fs_ports::PortKind::MechanicalForceVelocity,
                    false,
                ));
            }
        }
        let system = Self {
            models, candidates, exact, offsets, connections, columns,
            config, dt, old_q: vec![0.0; count], free_delta: vec![0.0; count],
            forces: vec![0.0; count],
            reactions: vec![0.0; links],
            staged_frame: frame.clone(),
            frame,
            staged_ports: ports.clone(), ports,
        };
        limit("total initial energy", system.total_energy_j()?, config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(system)
    }

    /// Accepted components in original order; callers cannot mutate their basis.
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] { &self.models }
    /// Flattened force layout: component order, then each component's mode order.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.forces.len() }
    /// Shared mechanical sample period [s].
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 { self.dt }
    /// Completed samples only. Failed trials consume no physical time.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.frame.sample }
    /// Last complete diagnostic, absent before the first step.
    #[must_use]
    pub fn last_frame(&self) -> Option<&CoupledModalFrame> {
        (self.frame.sample != 0).then_some(&self.frame)
    }
    /// Current oscillator and spring storage, including unexcited/preloaded input.
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let mut energy = component_energy(&self.models)?;
        for (column, link) in self.columns.iter().zip(&self.connections) {
            let x = extension(&self.models, column, link.rest_extension_m)?;
            energy = finite(energy + 0.5 * link.stiffness_n_m * x * x)?;
        }
        Ok(energy)
    }

    /// Advance one held external-force sample. All components and the energy
    /// report publish together after reaction, component and energy gates pass.
    pub fn step(&mut self, external: &[f64]) -> Result<&CoupledModalFrame, ModalCouplingError> {
        self.step_inner(external, None)
    }
    /// The same transaction with cancellation polls at component/link boundaries.
    /// A cancelled trial may be retried with unchanged accepted state and inputs.
    pub fn step_under_gate(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<&CoupledModalFrame, ModalCouplingError>
    {
        self.step_inner(external, Some(gate))
    }

    fn step_inner(&mut self, external: &[f64], gate: Option<&CancelGate>)
        -> Result<&CoupledModalFrame, ModalCouplingError>
    {
        self.stage_inner(external, gate)?;
        poll(gate)?;
        self.publish_staged();
        Ok(&self.frame)
    }

    // Contact extensions inspect the complete candidate before publication.
    // Scratch may change on refusal; accepted models and diagnostics never do.
    fn stage_inner(&mut self, external: &[f64], gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        exact::stage(self,external,gate)
    }

    fn publish_staged(&mut self) {
        std::mem::swap(&mut self.models, &mut self.candidates);
        std::mem::swap(&mut self.frame, &mut self.staged_frame);
        std::mem::swap(&mut self.ports, &mut self.staged_ports);
    }

    // Unrestrained full-network prediction is scratch. Contact/friction may
    // suppress it before any physical endpoint or force budget is applied.
    fn prepare_prediction(
        &mut self,
        external: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalCouplingError> {
        poll(gate)?;
        if external.len() != self.mode_count() || external.iter().any(|x| !x.is_finite()) {
            return Err(invalid("external forces must be finite and match all component modes"));
        }
        if let Some(exact) = &self.exact {
            let (record, x, _) = exact.advance(&self.models, external, self.dt)?;
            for k in 0..self.mode_count() {
                self.old_q[k] = x[2 * k];
                self.free_delta[k] = finite(record.x[2 * k] - x[2 * k])?;
            }
        } else {
            let mut k = 0;
            for model in &self.models {
                poll(gate)?;
                for (mode, state) in model.states().iter().enumerate() {
                    self.old_q[k] = state.displacement_m_sqrt_kg;
                    self.free_delta[k] = finite(
                        model
                            .advance_exact_zoh_step(mode, *state, external[k])
                            .displacement_m_sqrt_kg
                            - self.old_q[k],
                    )?;
                    k += 1;
                }
            }
        }
        self.forces.copy_from_slice(external);
        poll(gate)?;
        Ok(())
    }

    fn linear_response(
        &self,
        column: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ModalCouplingError> {
        poll(Some(gate))?;
        if let Some(exact) = &self.exact {
            return exact.response(column, self.dt);
        }
        if column.len() != self.mode_count() {
            return Err(invalid("response column must match every original mode"));
        }
        let mut result = Vec::with_capacity(column.len());
        let mut k = 0;
        for model in &self.models {
            poll(Some(gate))?;
            for mode in 0..model.modes().len() {
                result.push(finite(
                    column[k]
                        * model
                            .advance_exact_zoh_step(mode, ModalAcousticState::default(), 1.0)
                            .displacement_m_sqrt_kg,
                )?);
                k += 1;
            }
        }
        Ok(result)
    }
}

fn connection_layout(
    models: &[ModalAcousticTimeModel], links: usize,
    config: ModalCouplingConfig, gate: &CancelGate,
) -> Result<(Vec<usize>, usize, f64), ModalCouplingError> {
    poll(Some(gate))?;
    validate_config(config)?;
    if models.is_empty() || models.len() > config.max_modes || links > config.max_connections {
        return Err(invalid("component/connection count exceeds its nonempty admission"));
    }
    let dt = models[0].sample_period_s();
    let mut offsets = vec![0_usize];
    let mut count = 0_usize;
    for model in models {
        if model.sample_period_s().to_bits() != dt.to_bits() {
            return Err(invalid("every component must use the same sample clock"));
        }
        count = count.checked_add(model.modes().len()).ok_or_else(|| invalid("mode count overflow"))?;
        if count > config.max_modes { return Err(invalid("total modal count exceeds max_modes")); }
        offsets.push(count);
    }
    let terms = count.checked_mul(links + 1).and_then(|n| n.checked_mul(links + 1))
        .ok_or_else(|| invalid("connection setup work overflow"))?;
    if terms > config.max_setup_terms { return Err(invalid("connection setup work exceeds max_setup_terms")); }
    Ok((offsets, count, dt))
}

fn validate_config(c: ModalCouplingConfig) -> Result<(), ModalCouplingError> {
    if c.max_modes == 0 || c.max_modes > 4096 || c.max_connections > 64
        || !c.nyquist_guard_fraction.is_finite() || c.nyquist_guard_fraction <= 0.0 || c.nyquist_guard_fraction >= 1.0
        || [c.maximum_total_energy_j, c.maximum_abs_pressure_pa, c.maximum_abs_connection_force_n,
            c.solve_relative_tolerance].iter().any(|x| !x.is_finite() || *x <= 0.0)
        || c.solve_relative_tolerance >= 1.0 || !c.energy_absolute_tolerance_j.is_finite()
        || c.energy_absolute_tolerance_j < 0.0 || !c.energy_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.energy_relative_tolerance)
        || (c.energy_absolute_tolerance_j == 0.0 && c.energy_relative_tolerance == 0.0) {
        return Err(invalid("invalid modal coupling budgets, Nyquist guard or residual tolerances"));
    }
    Ok(())
}
fn connection_matrix(columns: &[Vec<f64>], compliance: &[f64], roots: &[f64], gate: Option<&CancelGate>)
    -> Result<Vec<f64>, ModalCouplingError>
{
    let n = columns.len();
    let mut a = vec![0.0; n*n];
    for row in 0..n {
        poll(gate)?;
        for col in 0..=row {
            let mut sum = 0.0;
            for ((&left, &right), &d) in columns[row].iter().zip(&columns[col]).zip(compliance) {
                sum = finite(sum + left * d * right)?;
            }
            let value = finite(roots[row] * sum * roots[col] + if row == col { 1.0 } else { 0.0 })?;
            a[row*n+col] = value;
            a[col*n+row] = value;
        }
    }
    Ok(a)
}
fn check_solve(a: &[f64], x: &[f64], b: &[f64], tolerance: f64) -> Result<f64, ModalCouplingError> {
    let mut worst = 0.0_f64;
    for row in 0..b.len() {
        let mut applied = 0.0;
        let mut scale = b[row].abs();
        for col in 0..b.len() {
            let term = finite(a[row*b.len()+col] * x[col])?;
            applied = finite(applied + term)?;
            scale = finite(scale + term.abs())?;
        }
        let residual = finite(b[row] - applied)?.abs();
        worst = worst.max(if scale == 0.0 { 0.0 } else { residual / scale });
    }
    if worst > tolerance { return Err(ModalCouplingError::SolveResidual { relative: worst, tolerance }); }
    Ok(worst)
}
fn component_energy(models: &[ModalAcousticTimeModel]) -> Result<f64, ModalCouplingError> {
    let mut sum = 0.0;
    for model in models {
        for (mode, state) in model.modes().iter().zip(model.states()) {
            let v = state.velocity_m_sqrt_kg_per_s;
            let q = state.displacement_m_sqrt_kg;
            let w = mode.angular_frequency_rad_s;
            sum = finite(sum + f64::midpoint(v*v, w*w*q*q))?;
        }
    }
    Ok(sum)
}
fn extension(models: &[ModalAcousticTimeModel], column: &[f64], rest: f64) -> Result<f64, ModalCouplingError> {
    let mut sum = 0.0;
    for (state, &shape) in models.iter().flat_map(|m| m.states()).zip(column) {
        sum = finite(sum + shape * state.displacement_m_sqrt_kg)?;
    }
    finite(sum - rest)
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, ModalCouplingError> {
    a.iter().zip(b).try_fold(0.0, |s, (a,b)| finite(s+a*b))
}
fn finite(value: f64) -> Result<f64, ModalCouplingError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("derived coupling value is not finite")) }
}
fn limit(what: &'static str, value: f64, limit: f64) -> Result<(), ModalCouplingError> {
    if value > limit { Err(ModalCouplingError::Budget { what, value, limit }) } else { Ok(()) }
}
fn poll(gate: Option<&CancelGate>) -> Result<(), ModalCouplingError> {
    if gate.is_some_and(CancelGate::is_requested) { Err(ModalCouplingError::Cancelled) } else { Ok(()) }
}
fn invalid(what: &'static str) -> ModalCouplingError { ModalCouplingError::Invalid { what } }

/// Static preload and local implicit equilibrium sensitivities.
pub mod equilibrium;
/// Scheduled rendering of connected modal components.
pub mod render;

/// A bilateral network with one two-body unilateral power-law contact.
pub mod contact;

crate::port_adoption::impl_ported_module!(CoupledModalSystem, ports);
