//! Two-way elastic/viscous connections between retained modal components.
//!
//! Component evolution stays in `modal_acoustic_time`'s exact held-force
//! transition. Only the connection reactions are solved here. For extension
//! x = B^T q - rest, a connection applies lambda = -k (x0+x1)/2
//! - c (x1-x0)/dt through the SAME signed shape column B. Its work is
//! -delta(k*x*x/2) - c*(x1-x0)^2/dt; a receiver is not a one-way audio filter.
//!
//! With the component displacement compliance D for this step, h=k/2+c/dt,
//! lambda=-sqrt(h)*y, the small symmetric system is
//! (I + sqrt(h) B^T D B sqrt(h)) y = k*x0/sqrt(h) + sqrt(h)*delta_x_free.
//! Its Cholesky factor is prepared once by fs-la, not rebuilt per sample.
//! This is a finite-step port coupling, NOT the exact continuous coupled-system
//! exponential. Time refinement is still necessary. Pressure transfers retain
//! their original narrow-band, read-only scope; radiation loading is not added.

use fs_exec::CancelGate;
use fs_la::factor::{Cholesky, cholesky};

use crate::modal_acoustic_time::{
    ModalAcousticState, ModalAcousticTimeError, ModalAcousticTimeModel,
    ModalAcousticWorkspace, advance_exact_zoh,
};

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
    /// Sum of the component stepper's work-minus-energy viscous losses [J].
    pub component_dissipation_j: f64,
    /// Nonnegative c*(delta_extension)^2/dt of the connections [J].
    pub connection_dissipation_j: f64,
    /// Combined storage change + dissipation - external work [J].
    pub energy_residual_j: f64,
    /// Actual absolute + relative energy tolerance used by the gate [J].
    pub energy_tolerance_j: f64,
    /// Independently recomputed connection equation residual.
    pub solve_relative_residual: f64,
    /// Signed force on the LEFT side of each connection [N].
    pub connection_forces_n: Vec<f64>,
}

/// A fixed network with two-way reactions and transactional accepted states.
/// Models, their individual budgets, damping and pressure transfers are retained.
/// All state and solve buffers are allocated at construction. This is not a
/// measured real-time or experimental-validity claim.
pub struct CoupledModalSystem {
    models: Vec<ModalAcousticTimeModel>,
    candidates: Vec<ModalAcousticTimeModel>,
    workspaces: Vec<ModalAcousticWorkspace>,
    offsets: Vec<usize>,
    connections: Vec<ModalConnection>,
    columns: Vec<Vec<f64>>,
    roots: Vec<f64>,
    spring_over_root: Vec<f64>,
    matrix: Vec<f64>,
    factor: Cholesky,
    config: ModalCouplingConfig,
    dt: f64,
    old_q: Vec<f64>,
    free_delta: Vec<f64>,
    forces: Vec<f64>,
    rhs: Vec<f64>,
    solution: Vec<f64>,
    reactions: Vec<f64>,
    frame: CoupledModalFrame,
    staged_frame: CoupledModalFrame,
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
        poll(Some(gate))?;
        validate_config(config)?;
        if models.is_empty() || models.len() > config.max_modes || connections.len() > config.max_connections {
            return Err(invalid("component/connection count exceeds its nonempty admission"));
        }
        let dt = models[0].sample_period_s();
        let mut offsets = vec![0_usize];
        let mut count = 0_usize;
        for model in &models {
            if model.sample_period_s().to_bits() != dt.to_bits() {
                return Err(invalid("every component must use the same sample clock"));
            }
            count = count.checked_add(model.modes().len()).ok_or_else(|| invalid("mode count overflow"))?;
            if count > config.max_modes { return Err(invalid("total modal count exceeds max_modes")); }
            offsets.push(count);
        }
        let links = connections.len();
        let terms = count.checked_mul(links + 1).and_then(|n| n.checked_mul(links + 1))
            .ok_or_else(|| invalid("connection setup work overflow"))?;
        if terms > config.max_setup_terms { return Err(invalid("connection setup work exceeds max_setup_terms")); }
        let mut columns = Vec::with_capacity(links);
        let mut roots = Vec::with_capacity(links);
        let mut spring_over_root = Vec::with_capacity(links);
        for link in &connections {
            poll(Some(gate))?;
            if !link.stiffness_n_m.is_finite() || link.stiffness_n_m < 0.0
                || !link.damping_n_s_m.is_finite() || link.damping_n_s_m < 0.0
                || !link.rest_extension_m.is_finite() {
                return Err(invalid("connection stiffness/damping must be finite nonnegative; rest must be finite"));
            }
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
            let h = finite(0.5 * link.stiffness_n_m + link.damping_n_s_m / dt)?;
            if h == 0.0 && (link.stiffness_n_m > 0.0 || link.damping_n_s_m > 0.0) {
                return Err(invalid("positive connection coefficient is not representable"));
            }
            let root = h.sqrt();
            spring_over_root.push(if root == 0.0 { 0.0 } else { finite(link.stiffness_n_m / root)? });
            roots.push(root);
            columns.push(column);
        }
        let mut compliance = Vec::with_capacity(count);
        let mut upper_squared = 0.0_f64;
        for model in &models {
            poll(Some(gate))?;
            for &mode in model.modes() {
                upper_squared = upper_squared.max(finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s)?);
                let response = advance_exact_zoh(mode, ModalAcousticState::default(), 1.0, dt).displacement_m_sqrt_kg;
                if !response.is_finite() || response < 0.0 {
                    return Err(invalid("modal held-force displacement response must be finite nonnegative"));
                }
                compliance.push(response);
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
        let matrix = connection_matrix(&columns, &compliance, &roots, Some(gate))?;
        let factor = cholesky(&matrix, links).map_err(ModalCouplingError::Factor)?;
        let candidates = models.clone();
        let workspaces = models.iter().map(ModalAcousticWorkspace::new).collect();
        let frame = CoupledModalFrame {
            sample: 0, observer_pressure_pa: 0.0, modal_energy_j: 0.0, connection_energy_j: 0.0,
            external_work_j: 0.0, component_dissipation_j: 0.0, connection_dissipation_j: 0.0,
            energy_residual_j: 0.0, energy_tolerance_j: 0.0, solve_relative_residual: 0.0,
            connection_forces_n: vec![0.0; links],
        };
        let system = Self {
            models, candidates, workspaces, offsets, connections, columns, roots, spring_over_root,
            matrix, factor, config, dt, old_q: vec![0.0; count], free_delta: vec![0.0; count],
            forces: vec![0.0; count], rhs: vec![0.0; links], solution: vec![0.0; links],
            reactions: vec![0.0; links],
            staged_frame: frame.clone(),
            frame,
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
        poll(gate)?;
        if external.len() != self.mode_count() || external.iter().any(|x| !x.is_finite()) {
            return Err(invalid("external forces must be finite and match all component modes"));
        }
        let sample = self.frame.sample.checked_add(1).ok_or_else(|| invalid("coupled sample clock overflow"))?;
        let before = self.total_energy_j()?;
        let relative = self.prepare_forces(external, gate, true)?;
        let mut pressure = 0.0;
        let mut modal_after = 0.0;
        let mut component_loss = 0.0;
        let mut work = 0.0;
        for i in 0..self.models.len() {
            poll(gate)?;
            self.candidates[i].restore_states(self.models[i].states())
                .map_err(|source| ModalCouplingError::Component { component: i, source })?;
            let range = self.offsets[i]..self.offsets[i+1];
            let frame = self.candidates[i].step_into(&self.forces[range.clone()], &mut self.workspaces[i])
                .map_err(|source| ModalCouplingError::Component { component: i, source })?;
            pressure = finite(pressure + frame.observer_pressure_pa)?;
            modal_after = finite(modal_after + frame.total_modal_energy_j)?;
            component_loss = finite(component_loss + frame.viscous_dissipation_j)?;
            for (index, state) in range.zip(self.candidates[i].states()) {
                let delta = finite(state.displacement_m_sqrt_kg - self.old_q[index])?;
                self.free_delta[index] = delta;
                work = finite(work + external[index] * delta)?;
            }
        }
        let mut spring_after = 0.0;
        let mut connection_loss = 0.0;
        for (column, link) in self.columns.iter().zip(&self.connections) {
            poll(gate)?;
            let x = extension(&self.candidates, column, link.rest_extension_m)?;
            spring_after = finite(spring_after + 0.5 * link.stiffness_n_m * x * x)?;
            let delta = dot(column, &self.free_delta)?;
            connection_loss = finite(connection_loss + link.damping_n_s_m * (delta / self.dt) * delta)?;
        }
        let after = finite(modal_after + spring_after)?;
        limit("total energy", after, self.config.maximum_total_energy_j)?;
        limit("observer pressure", pressure.abs(), self.config.maximum_abs_pressure_pa)?;
        let residual = finite((after - before) + component_loss + connection_loss - work)?;
        let scale = before.max(after).max(work.abs()).max(component_loss.abs()).max(connection_loss);
        let tolerance = finite(self.config.energy_absolute_tolerance_j + self.config.energy_relative_tolerance * scale)?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance });
        }
        poll(gate)?;
        self.staged_frame.sample = sample;
        self.staged_frame.observer_pressure_pa = pressure;
        self.staged_frame.modal_energy_j = modal_after;
        self.staged_frame.connection_energy_j = spring_after;
        self.staged_frame.external_work_j = work;
        self.staged_frame.component_dissipation_j = component_loss;
        self.staged_frame.connection_dissipation_j = connection_loss;
        self.staged_frame.energy_residual_j = residual;
        self.staged_frame.energy_tolerance_j = tolerance;
        self.staged_frame.solve_relative_residual = relative;
        self.staged_frame.connection_forces_n.copy_from_slice(&self.reactions);
        Ok(())
    }

    fn publish_staged(&mut self) {
        std::mem::swap(&mut self.models, &mut self.candidates);
        std::mem::swap(&mut self.frame, &mut self.staged_frame);
    }

    // The single bilateral-reaction owner, shared by accepted samples and the
    // contact-free affine prediction. Only FINAL reactions face physical caps:
    // contact can suppress an otherwise over-budget unrestrained prediction.
    fn prepare_forces(&mut self, external: &[f64], gate: Option<&CancelGate>, enforce_limits: bool)
        -> Result<f64, ModalCouplingError>
    {
        poll(gate)?;
        if external.len() != self.mode_count() || external.iter().any(|x| !x.is_finite()) {
            return Err(invalid("external forces must be finite and match all component modes"));
        }
        let mut at = 0;
        for model in &self.models {
            poll(gate)?;
            for (&mode, &state) in model.modes().iter().zip(model.states()) {
                self.old_q[at] = state.displacement_m_sqrt_kg;
                // Unconstrained prediction is scratch, not an accepted component
                // step: a stiff connection may suppress a large free prediction.
                self.free_delta[at] = finite(advance_exact_zoh(mode, state, external[at], self.dt)
                    .displacement_m_sqrt_kg - state.displacement_m_sqrt_kg)?;
                at += 1;
            }
        }
        for (j, (column, link)) in self.columns.iter().zip(&self.connections).enumerate() {
            poll(gate)?;
            let x = finite(dot(column, &self.old_q)? - link.rest_extension_m)?;
            self.rhs[j] = finite(self.spring_over_root[j] * x + self.roots[j] * dot(column, &self.free_delta)?)?;
        }
        self.solution.copy_from_slice(&self.rhs);
        self.factor.solve(&mut self.solution);
        let relative = check_solve(&self.matrix, &self.solution, &self.rhs, self.config.solve_relative_tolerance)?;
        self.forces.copy_from_slice(external);
        for j in 0..self.connections.len() {
            poll(gate)?;
            let reaction = finite(-self.roots[j] * self.solution[j])?;
            if enforce_limits {
                limit("connection force", reaction.abs(), self.config.maximum_abs_connection_force_n)?;
            }
            self.reactions[j] = reaction;
            // Zero reactions leave the original external force bits alone.
            if reaction != 0.0 {
                for (force, &b) in self.forces.iter_mut().zip(&self.columns[j]) {
                    *force = finite(*force + b * reaction)?;
                }
            }
        }
        Ok(relative)
    }

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

mod equilibrium;
/// Scheduled rendering of connected modal components.
pub mod render;

/// A bilateral network with one two-body unilateral power-law contact.
pub mod contact;
