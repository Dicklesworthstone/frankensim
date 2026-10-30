//! Implicit driven dry friction on the existing two-way modal network.
//!
//! The attachment coordinate is x=B^T q (left minus optional right). A held
//! traction F acts through B; the prescribed surface advances by u*dt. The
//! SAME displacement increment defines slip and work: d=u*dt-B^T delta_q.
//! Sticking solves d=0 inside the static cone. Sliding solves F=mu(|d|/dt) N
//! with the sign of d, using fs-tribo's kinetic law, not a velocity regularizer.
//! The network supplies its complete bilateral held-force response, so every
//! connected body reacts back during this solve. Its original time stepper,
//! component/connection limits and state transaction remain in charge.
//!
//! This is step-averaged sticking, not an exact endpoint velocity constraint.
//! One tangential coordinate and a prescribed nonnegative normal load are
//! supported. Where a Stribeck or hyperbolic law meets the load line more than
//! once, the contact keeps the branch it occupies (McIntyre, Schumacher &
//! Woodhouse 1983, through `stribeck_friction::StribeckLoadLine`), and every
//! frame discloses the root count and the committed root. Other laws keep the
//! scalar bisection. No normal collision, finite contact patch, thermal/wear
//! evolution, experimental validation or hard-real-time claim is implied. Input law coefficients retain caller authority. A nonzero drive
//! speed is an explicit work source; zero drive gives a passive internal brake.

use fs_exec::CancelGate;
use fs_tribo::{FrictionLaw, TriboError};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use crate::stribeck_friction::{ContactBranch, ImplicitFrictionError, StribeckLoadLine};
use super::effective_compliance;
use super::super::{CoupledModalSystem, ModalAttachment, ModalCouplingError,
    dot, finite, invalid, limit, poll};

/// Signed tangential attachment; a missing right side is a prescribed surface.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalFrictionPort {
    /// Positive side of x, receiving +B_left F.
    pub left: ModalAttachment,
    /// Optional negative side, receiving -B_right F in the SAME solve.
    pub right: Option<ModalAttachment>,
    /// Unmodified caller-declared fs-tribo law, without inferred source authority.
    pub law: FrictionLaw,
}

/// External controls held over one accepted mechanical sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionDrive {
    /// Signed prescribed relative surface speed [m/s]. Zero is an undriven brake.
    pub speed_m_s: f64,
    /// Prescribed compressive load [N]. For a bound adhesive patch this must
    /// equal its current repulsive load held over this sample; zero need not
    /// release junction shear.
    pub normal_force_n: f64,
}

/// Caller-selected physical limits and nonlinear-work budget; no clipping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalFrictionConfig {
    /// At most this many bisections, in 1..=128.
    pub max_iterations: usize,
    /// Maximum absolute accepted tangential traction [N].
    pub maximum_force_n: f64,
    /// Maximum absolute accepted step-averaged slip speed [m/s].
    pub maximum_slip_speed_m_s: f64,
    /// Absolute force-law residual tolerance [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative force-law residual tolerance, in (0,1).
    pub force_relative_tolerance: f64,
    /// Absolute displacement mismatch allowed on the sticking branch [m].
    pub sticking_tolerance_m: f64,
}

/// Resolved discrete branch, never inferred from a small velocity threshold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModalFrictionRegime {
    /// The dry normal load or adhesive real area was zero; the network still evolves.
    Released,
    /// The displacement constraint is satisfied within its declared tolerance.
    Sticking,
    /// The fs-tribo kinetic law is satisfied at actual step-averaged slip.
    Sliding,
    /// An explicitly bound material layer carries reversible elastic shear.
    /// Finite elastic relative motion is permitted even without plastic flow.
    Shearing {
        /// Whether this interval produced irreversible plastic slip.
        plastic_flow: bool,
    },
}

/// Network, constitutive and step-averaged sticking failures stay distinct.
#[derive(Debug)]
pub enum ModalFrictionError {
    /// An original numerical owner or whole-system budget refused the trial.
    Network(ModalCouplingError),
    /// The unchanged fs-tribo law refused its parameters or slip speed.
    Law(TriboError),
    /// Original adhesive source identity or geometry admission refused.
    AdhesiveSource(fs_contact::interface_binding::adhesive_shear::AdhesiveShearError),
    /// Actual endpoint displacements failed the sticking constraint.
    Sticking { /// Signed mismatch [m].
        residual_m: f64, /// Explicit tolerance [m].
        tolerance_m: f64 },
    /// A bounded static constitutive solve did not close the actual force law.
    NonlinearEquilibrium {
        /// Force mismatch at the worst normalized residual [N/sqrt(kg)].
        residual: f64,
        /// Its absolute-plus-relative allowance [N/sqrt(kg)].
        tolerance: f64,
        /// Constitutive trial evaluations spent, including line-search trials.
        iterations: usize,
    },
    /// The active-set rescue (frankensim-o7qbp) could not resolve a consistent
    /// joint stick/slip partition within its pivot budget. Carries the pivots
    /// spent against the budget, the last active set tried (S/+/- per port) and
    /// the worst Karush-Kuhn-Tucker violation (port, condition, size).
    ActiveSetExhausted {
        /// Pivots the active-set walk spent.
        pivots: usize,
        /// Pivot budget it was allowed.
        budget: usize,
        /// Last active set tried, one char per port (S stuck, + slide+, - slide-).
        last_active_set: String,
        /// Port of the worst remaining KKT violation.
        worst_port: usize,
        /// Which friction condition that port broke.
        worst_condition: &'static str,
        /// Magnitude of that violation.
        worst_size: f64,
    },
    /// The active-set rescue was not attempted because an obstacle or shear layer
    /// is present on a port; the joint solve refused on its original terms.
    ActiveSetDeclined {
        /// The joint sweep's own failure at its last sweep: what did not converge.
        cause: Box<ModalFrictionError>,
    },
}
impl From<ModalCouplingError> for ModalFrictionError {
    fn from(value: ModalCouplingError) -> Self { Self::Network(value) }
}
impl From<TriboError> for ModalFrictionError {
    fn from(value: TriboError) -> Self { Self::Law(value) }
}
impl From<fs_contact::interface_binding::adhesive_shear::AdhesiveShearError> for ModalFrictionError {
    fn from(value: fs_contact::interface_binding::adhesive_shear::AdhesiveShearError) -> Self { Self::AdhesiveSource(value) }
}
impl core::fmt::Display for ModalFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Network(e) => write!(f, "friction network: {e}"),
            Self::Law(e) => write!(f, "friction law: {e}"),
            Self::AdhesiveSource(e) => write!(f, "friction source: {e}"),
            Self::Sticking { residual_m, tolerance_m } => write!(f,
                "sticking mismatch {residual_m:e} m exceeds {tolerance_m:e} m"),
            Self::NonlinearEquilibrium { residual, tolerance, iterations } => write!(f,
                "nonlinear static force residual {residual:e} exceeds {tolerance:e} N/sqrt(kg) after {iterations} trials"),
            Self::ActiveSetExhausted { pivots, budget, last_active_set, worst_port, worst_condition, worst_size } => write!(f,
                "active-set rescue exhausted: {pivots}/{budget} pivots, last active set {last_active_set}, \
                 worst KKT violation at port {worst_port} ({worst_condition}) size {worst_size:e}"),
            Self::ActiveSetDeclined { cause } => write!(f,
                "active-set rescue not attempted: an obstacle or shear layer is present on a port; \
                 the joint sweep's last failure: {cause}"),
        }
    }
}
impl std::error::Error for ModalFrictionError {}

/// A complete accepted energy/work window, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionModalFrame {
    /// Original network's one-based accepted sample clock.
    pub sample: u64,
    /// Original read-only pressure sum, or the explicit acceleration observation [Pa].
    pub observer_pressure_pa: f64,
    /// Modal, bilateral-spring, material-memory and bound nonlinear storage [J].
    pub network_energy_j: f64,
    /// Work of external modal forces, excluding friction [J].
    pub external_work_j: f64,
    /// Work of the prescribed moving surface, F*u*dt [J]. May be negative.
    pub drive_work_j: f64,
    /// Original component, bilateral-dashpot and material-relaxation loss [J].
    pub network_dissipation_j: f64,
    /// Nonnegative F*d on sliding; zero for sticking/released [J].
    pub friction_dissipation_j: f64,
    /// Storage change + losses - external and surface work [J].
    pub energy_residual_j: f64,
    /// Original network's scaled energy tolerance [J].
    pub energy_tolerance_j: f64,
    /// Signed traction on the positive attachment [N].
    pub traction_n: f64,
    /// Actual u*dt - B^T delta_q [m], including sticking roundoff.
    pub slip_distance_m: f64,
    /// Actual signed kinetic residual [N]; zero on other branches.
    pub constitutive_residual_n: f64,
    /// Resolved branch for this step.
    pub regime: ModalFrictionRegime,
    /// Bisection refinements used; zero for a closed-form branch or endpoint.
    /// The branch-continuation owner bisects each root to adjacent floats and
    /// refuses when that exceeds `max_iterations`.
    pub iterations: usize,
    /// Root disclosure of a Stribeck or hyperbolic step; `None` for laws
    /// without a disclosed root rule and for a released contact.
    pub roots: Option<FrictionRoots>,
}

/// Intersections of the step's load line with the friction graph (the static
/// cone included) and the zero-based index of the committed one, ordered by
/// slip speed with the cone first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrictionRoots {
    /// Intersections found, 1..=4; four requires a static/kinetic junction drop.
    pub found: u8,
    /// Index of the committed intersection.
    pub selected: u8,
}

/// Retained physical network with one implicitly solved tangential port.
/// Admission adds no stored energy and can preserve an already vibrating network.
/// Refusal/cancellation publishes neither candidate state nor a new frame.
pub struct FrictionModalSystem {
    network: CoupledModalSystem,
    port: ModalFrictionPort,
    config: ModalFrictionConfig,
    column: Vec<f64>,
    compliance: f64,
    forces: Vec<f64>,
    last: Option<FrictionModalFrame>,
    acceleration_pressure: Option<Vec<f64>>,
    interface: Option<interface::DryFrictionBinding>,
    adhesive: Option<adhesive::AdhesiveModalBinding>,
    nonlinear: Option<nonlinear::NonlinearFrictionState>,
    // McIntyre-Schumacher-Woodhouse branch memory: committed, and staged by
    // the trial that may still be refused or repeated.
    branch: ContactBranch,
    staged_branch: ContactBranch,
    load_line: Option<StribeckLoadLine>,
}

impl FrictionModalSystem {
    /// Condense the admitted bilateral response once; never substitute a fixed
    /// receiver or an uncoupled effective mass for the real mechanical network.
    pub fn new(network: CoupledModalSystem, port: ModalFrictionPort,
        config: ModalFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        let column = scalar::admit_port(&network, &port, config)?;
        let compliance = construction_compliance(&network, &column, gate)?;
        poll(Some(gate))?;
        Ok(Self { forces: vec![0.0; network.mode_count()], network, port,
            config, column, compliance, last: None, acceleration_pressure: None, interface: None, adhesive: None, nonlinear: None,
            branch: ContactBranch::Fresh, staged_branch: ContactBranch::Fresh, load_line: None })
    }

    /// Explicitly prepare the SAME full linear operator before playback.
    /// Refresh its exact held-force compliance; preserve the original friction
    /// law, material interface, pressure weights and component initial states.
    pub fn with_prepared_linear_propagation(
        mut self, budget: super::super::PreparedModalBudget, gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.network = self.network.with_prepared_linear_propagation(budget, gate)?;
        self.compliance = construction_compliance(&self.network, &self.column, gate)?;
        poll(Some(gate))?;
        Ok(self)
    }

    /// Bind constant-loss dashpots to the prepared coupled owner before playback
    /// ([`CoupledModalSystem::with_constant_loss_connections`]) and refresh its
    /// exact held-force compliance. The friction law, material interface,
    /// pressure weights and initial states are unchanged.
    pub fn with_constant_loss_connections(
        mut self, constant_loss: &[super::super::ConstantLossConnection], gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.network = self.network.with_constant_loss_connections(constant_loss, gate)?;
        self.compliance = construction_compliance(&self.network, &self.column, gate)?;
        poll(Some(gate))?;
        Ok(self)
    }

    /// The coupled owner's constant-loss connections, in connection order.
    pub fn constant_loss_connections(&self) -> Vec<super::super::ConstantLossConnection> {
        self.network.constant_loss_connections()
    }

    /// Bind genuine constitutive memory to the unchanged coupled owner before
    /// playback. Refresh the SAME friction force response, including material
    /// feedback; existing interfaces, nonlinear storage and pressure weights
    /// remain attached. See [`CoupledModalSystem::with_material_relaxation`].
    pub fn with_material_relaxation(
        mut self, branches: Vec<fs_phs::RelaxationBranch>, initial_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.network = self.network.with_material_relaxation(branches, initial_memory, gate)?;
        self.compliance = construction_compliance(&self.network, &self.column, gate)?;
        limit("complete initial material/nonlinear energy", self.total_energy_j()?,
            self.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(self)
    }

    /// Compose additional material participants in the same force response.
    /// Preserve existing history, interfaces, nonlinear storage and observers;
    /// refresh the bow compliance after the complete exact operator changes.
    pub fn with_additional_material_relaxation(
        mut self, branches: Vec<fs_phs::RelaxationBranch>, initial_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.network = self.network.with_additional_material_relaxation(branches, initial_memory, gate)?;
        self.compliance = construction_compliance(&self.network, &self.column, gate)?;
        limit("complete initial material/nonlinear energy", self.total_energy_j()?,
            self.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(self)
    }

    /// Unmodified material projections and coefficients, in memory order.
    pub fn with_plane_stress_material_points(mut self,
        points: Vec<super::super::ModalPlaneStressPoint>, gate:&CancelGate,
    ) -> Result<Self,ModalFrictionError> {
        self.network=self.network.with_plane_stress_material_points(points,gate)?;
        Ok(self)
    }

    /// Unmodified material projections and coefficients, in memory order.
    #[must_use]
    pub fn material_relaxation_branches(&self) -> &[fs_phs::RelaxationBranch] {
        self.network.material_relaxation_branches()
    }
    /// Accepted energy-normalized internal strains, including prior unloaded history.
    #[must_use]
    pub fn material_memory(&self) -> &[f64] { self.network.material_memory() }
    /// Current independent constitutive storage [J].
    pub fn material_energy_j(&self) -> Result<f64, ModalCouplingError> {
        self.network.material_energy_j()
    }
    /// Attach the SAME passive exterior load to the original held-force network
    /// and refresh the bow compliance including all acoustic reactions. Existing
    /// nonlinear potentials, source interfaces and pressure maps are retained.
    pub fn with_radiation_load(mut self,
        model: &crate::render::plate::impact::radiation::Model,
        source_modes: &[usize], maximum_poles: usize, gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.network = self.network.with_radiation_load(model, source_modes, maximum_poles, gate)?;
        self.compliance = construction_compliance(&self.network, &self.column, gate)?;
        poll(Some(gate))?;
        Ok(self)
    }
    /// Acoustic history participates in the same accepted sample transaction.
    #[must_use]
    pub fn radiation_state(&self) -> &[f64] { self.network.radiation_state() }
    /// Complete retained exterior impedance in the original physical port basis.
    #[must_use]
    pub fn radiation_model(&self) -> Option<&crate::render::plate::impact::radiation::Model> {
        self.network.radiation_model()
    }
    pub fn radiation_energy_j(&self) -> Result<f64, ModalCouplingError> {
        self.network.radiation_energy_j()
    }
    /// Last accepted linear-network account, with independent constitutive loss.
    #[must_use]
    pub fn coupled_frame(&self) -> Option<&super::super::CoupledModalFrame> {
        self.network.last_frame()
    }

    /// Explicit complete-operator preparation choice, if supplied.
    #[must_use]
    pub const fn prepared_budget(&self) -> Option<super::super::PreparedModalBudget> {
        self.network.prepared_budget()
    }

    /// Explicit construction purpose; harmonic-only networks retain this
    /// identity through material and radiation attachments.
    pub fn preparation(&self) -> super::super::ModalPreparation {
        self.network.preparation()
    }

    pub const fn is_harmonic_only(&self) -> bool { self.network.is_harmonic_only() }

    /// Accepted components, with their existing basis, budgets and state.
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] { self.network.components() }
    /// Original complete linear-network stiffness bound, including all springs.
    /// Independent friction/nonlinear trajectories need their own refinement.
    #[must_use]
    pub const fn linear_frequency_upper_bound_rad_s(&self) -> f64 {
        self.network.linear_frequency_upper_bound_rad_s()
    }
    /// Flattened external-force dimension, component then mode order.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.network.mode_count() }
    /// Original time step [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.network.sample_period_s() }
    /// Original accepted sample clock; refused trials never advance it.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.network.samples_rendered() }
    /// Current complete storage, including any explicitly bound nonlinear law [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        let linear = self.network.total_energy_j()?;
        if self.nonlinear.is_some() { finite(linear + self.nonlinear_energy_j()?) }
        else { Ok(linear) }
    }
    /// Settle external loads against the COMPLETE connected instrument with all
    /// friction contacts released. Uses the original static inverse stiffness;
    /// compression/bridge springs and any supported free coordinates remain
    /// coupled. Requires sample zero and zero initial Q/V, so existing vibration
    /// cannot be erased. No clock advances and no state publishes on refusal.
    /// Bound conservative storage participates through its exact gradient and
    /// the existing nonlinear force/iteration limits. Returned stored energy is
    /// an initial condition, NOT measured settling work.
    pub fn initialize_released_equilibrium(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalFrictionError>
    {
        if let Some(state) = &mut self.nonlinear {
            let report = state.initialize_equilibrium(&mut self.network, external, None, gate)?;
            let energy = report.total_energy_j;
            state.preload = Some(report);
            return Ok(energy);
        }
        self.network.initialize_static_equilibrium(external, gate).map_err(Into::into)
    }

    /// Last fully accepted friction-inclusive frame.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&FrictionModalFrame> { self.last.as_ref() }
    /// Original attachment and dry pressure-law component. When adhesion is
    /// bound, `adhesive_binding` retains its additional junction law and sources.
    #[must_use]
    pub const fn port(&self) -> &ModalFrictionPort { &self.port }

    /// Advance one held-input physical sample without a cancellation request.
    pub fn step(&mut self, external: &[f64], drive: FrictionDrive)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.step_inner(external, drive, None) }
    /// Same transaction with bounded cancellation polls at prediction/root/stage boundaries.
    pub fn step_under_gate(&mut self, external: &[f64], drive: FrictionDrive, gate: &CancelGate)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.step_inner(external, drive, Some(gate)) }

    fn step_inner(&mut self, external: &[f64], drive: FrictionDrive, gate: Option<&CancelGate>)
        -> Result<&FrictionModalFrame, ModalFrictionError>
    {
        let (accepted, nonlinear_frame) = self.stage_candidate(external, drive, gate)?;
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        Ok(self.publish_candidate(accepted, nonlinear_frame))
    }

    // Thermal coupling repeats this same physical proposal at endpoint
    // temperatures. Accepted mechanics and constitutive history stay untouched.
    fn stage_candidate(&mut self, external: &[f64], drive: FrictionDrive, gate: Option<&CancelGate>)
        -> Result<(FrictionModalFrame, Option<nonlinear::NonlinearFrictionFrame>), ModalFrictionError>
    {
        self.network.require_time_domain()?;
        let (accepted, nonlinear_frame) = if let Some(mut state) = self.nonlinear.take() {
            let trial = state.stage(self, external, drive, gate);
            self.nonlinear = Some(state);
            let frame = trial?;
            (frame.friction, Some(frame))
        } else { (self.stage_inner(external, drive, gate)?, None) };
        Ok((accepted, nonlinear_frame))
    }

    fn publish_candidate(&mut self, accepted: FrictionModalFrame,
        nonlinear_frame: Option<nonlinear::NonlinearFrictionFrame>) -> &FrictionModalFrame
    {
        self.network.staged_frame.observer_pressure_pa = accepted.observer_pressure_pa;
        self.network.publish_staged();
        self.branch = self.staged_branch;
        if let (Some(state), Some(frame)) = (&mut self.nonlinear, nonlinear_frame) {
            state.publish(frame);
        }
        self.last.insert(accepted)
    }

    // The nonlinear storage adapter may repeat this SAME held-force/friction
    // solve. Only the final jointly admitted candidate is ever published.
    fn stage_inner(&mut self, external: &[f64], drive: FrictionDrive, gate: Option<&CancelGate>)
        -> Result<FrictionModalFrame, ModalFrictionError>
    {
        poll(gate)?;
        if !drive.speed_m_s.is_finite() || !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0 {
            return Err(invalid("friction drive requires finite signed speed and nonnegative normal load").into());
        }
        interface::check_endpoint(self.interface.as_ref(), &self.column, &self.network.models, drive)?;
        let before = self.network.total_energy_j()?;
        self.network.prepare_prediction(external, gate)?;
        let free_increment=dot(&self.column,&self.network.free_delta)?;
        let drive_distance = finite(drive.speed_m_s*self.network.dt)?;
        let free_slip = finite(drive_distance-free_increment)?;
        let (traction, regime, iterations, roots) = self.solve_contact(free_slip, drive.normal_force_n, gate)?;
        self.forces.copy_from_slice(external);
        if traction != 0.0 {
            for (force, b) in self.forces.iter_mut().zip(&self.column) { *force = finite(*force+b*traction)?; }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let slip = finite(drive_distance-dot(&self.column, &self.network.free_delta)?)?;
        let (loss, constitutive_residual_n) = if let Some(binding) = &self.adhesive {
            binding.check(self.config, self.network.dt, traction, slip, regime, iterations)?
        } else { scalar::check(
            &self.port.law, self.config, self.network.dt, drive, traction, slip, regime, iterations,
        )? };
        interface::check_average(self.interface.as_ref(), drive, slip, self.network.dt)?;
        interface::check_endpoint(self.interface.as_ref(), &self.column, &self.network.candidates, drive)?;
        let observer_pressure_pa = self.staged_observer_pressure(gate)?;
        let frame = &self.network.staged_frame;
        let after = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let drive_work = finite(traction*drive_distance)?;
        let residual = finite((after-before)+network_loss+loss-external_work-drive_work)?;
        let scale = before.max(after).max(network_loss).max(loss).max(external_work.abs()).max(drive_work.abs());
        let tolerance = finite(self.network.config.energy_absolute_tolerance_j
            +self.network.config.energy_relative_tolerance*scale)?;
        if loss < 0.0 || residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance }.into());
        }
        let accepted = FrictionModalFrame {
            sample: frame.sample, observer_pressure_pa,
            network_energy_j: after, external_work_j: external_work, drive_work_j: drive_work,
            network_dissipation_j: network_loss, friction_dissipation_j: loss,
            energy_residual_j: residual, energy_tolerance_j: tolerance, traction_n: traction,
            slip_distance_m: slip, constitutive_residual_n, regime, iterations, roots,
        };
        poll(gate)?;
        Ok(accepted)
    }

    // Stribeck and hyperbolic laws take the branch-continuation owner on the
    // same held-force load line, v = v_free - (c/dt) F with v_free = free_slip/dt,
    // starting from the branch the contact occupies. Every other law, and a
    // released contact, keeps the scalar bisection unchanged.
    fn solve_contact(&mut self, free_slip: f64, normal: f64, gate: Option<&CancelGate>)
        -> Result<(f64, ModalFrictionRegime, usize, Option<FrictionRoots>), ModalFrictionError>
    {
        let continued = matches!(self.port.law, FrictionLaw::Stribeck { .. } | FrictionLaw::Hyperbolic { .. });
        if let Some(binding) = &self.adhesive {
            let contact = binding.tangential().contact();
            if normal != contact.repulsive_normal_force_n() {
                return Err(invalid("held adhesive input must match the current repulsive normal load").into());
            }
            if !continued || contact.area_m2() == 0.0 {
                let (traction, regime, iterations) = binding.solve(self.config, self.compliance,
                    self.network.dt, free_slip, normal, gate)?;
                self.staged_branch = ContactBranch::Fresh;
                return Ok((traction, regime, iterations, None));
            }
        } else if !continued || normal == 0.0 {
            self.staged_branch = ContactBranch::Fresh;
            let (traction, regime, iterations) = scalar::solve(&self.port.law, self.config, self.compliance, self.network.dt, free_slip, normal, gate)?;
            return Ok((traction, regime, iterations, None));
        }
        poll(gate)?;
        let mobility = finite(self.compliance/self.network.dt)?;
        let refusal = |error: ImplicitFrictionError<TriboError>| -> ModalFrictionError { match error {
            ImplicitFrictionError::Constitutive(cause) => cause.into(),
            ImplicitFrictionError::InvalidInput(what) => invalid(what).into(),
            ImplicitFrictionError::NonDissipative => invalid("set-valued friction would drive relative slip").into(),
            ImplicitFrictionError::NoConvergence => invalid("set-valued friction root did not resolve within max_iterations").into(),
        }};
        if !self.load_line.as_ref().is_some_and(|line| line.matches(normal, mobility)) {
            self.load_line = Some(if let Some(binding) = &self.adhesive {
                StribeckLoadLine::for_adhesive(binding.tangential(), mobility).map_err(refusal)?
            } else {
                StribeckLoadLine::new(&self.port.law, normal, mobility).map_err(refusal)?
            });
        }
        let line = self.load_line.as_ref().ok_or_else(|| invalid("missing set-valued load line"))?;
        let (step, iterations) = line.step_within(finite(free_slip/self.network.dt)?, self.branch,
            self.config.max_iterations).map_err(refusal)?;
        limit("friction traction", step.traction_n.abs(), self.config.maximum_force_n)?;
        self.staged_branch = step.branch;
        let regime = if step.branch == ContactBranch::Stuck { ModalFrictionRegime::Sticking } else { ModalFrictionRegime::Sliding };
        Ok((step.traction_n, regime, iterations, Some(FrictionRoots { found: step.roots, selected: step.selected })))
    }
}

// Analysis retains the original port/law but never consumes its held-step
// compliance. Every prediction and advancing path checks construction purpose
// before the private zero placeholder can participate in a contact solve.
fn construction_compliance(
    network: &CoupledModalSystem, column: &[f64], gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    poll(Some(gate))?;
    if network.is_harmonic_only() { return Ok(0.0); }
    effective_compliance(network, column, gate)
}

// Optional physical acceleration readout; it does not replace mechanical evolution.
mod observation;

/// Released-contact response of the original linear mechanical/acoustic owner.
pub mod harmonic;

mod scalar;

#[cfg(test)]
mod continuation_tests;

/// Several driven friction ports on the same retained mechanical network.
pub mod multiple;


/// Original dry-interface identities and domains on actual coupled contact motion.
pub mod interface;

/// Additive nonlinear storage on the original exact linear/friction network.
pub mod nonlinear;

/// Same-step friction, common-interface temperature and finite-body thermal evolution.
pub mod heat;

/// Adhesive junction shear with normal-owner geometry held over each sample.
pub mod adhesive;

/// Stateful tangential laws with reversible contact storage on the same network.
pub mod tangential;
