//! Multirate stepping with conserved interface accumulators
//! (bead `frankensim-rc-1109-root-3gc12.5.9`).
//!
//! Generalizes [`fs_phs::multirate::MultiratePair`] from two port-Hamiltonian
//! owners to arbitrary sealed [`System`] graphs across `G >= 1` rate groups:
//!
//! 1. [`System::plan_multirate`] assigns components to rate groups from each
//!    component's `limit_dt_s` and optional stiffness/substep overrides, keeps
//!    load-line loops intact in one rate group, builds per-group local Tarjan
//!    sub-blocks, and identifies all cross-rate interface outputs.
//! 2. [`System::step_multirate`] subcycles each rate group on its own local
//!    clock `dt_g = H / substeps_g`. When every group has `substeps == 1`, it
//!    dispatches directly to [`System::step`] for bitwise single-rate identity.
//! 3. Fast-to-slow interface exchange integrates every cross-rate output, port
//!    energy, entropy, stream mass, and carried inventory in double-double
//!    (`fs_math::dd::Dd`) over the macro window so nothing is created or lost
//!    at the rate boundary. Cyclic cross-rate stages close the macro window-
//!    average interface equations with safeguarded dense Newton (`fs_la::factor::lu`).
//! 4. Slow-to-fast coupling supports [`SlowInterpolation::ZeroOrderHold`] (exact
//!    work/heat/mass cancellation to `Dd` roundoff, matching `MultiratePair`)
//!    and [`SlowInterpolation::Linear`] (linear ramp across fast substeps whose
//!    exact quadrature energy discrepancy is explicitly computed in `Dd` and
//!    booked as numerical loss, never silently swallowed).
//! 5. [`System::step_multirate_to_events`] and mid-step component/boundary
//!    event handling synchronize every rate group to the exact same `Sum` clock
//!    with atomic checkpoint/rollback on any refusal or cancellation.

use core::any::Any;
use core::fmt;
use std::collections::BTreeMap;
use std::rc::Rc;

use fs_exec::Cx;
use fs_la::factor::lu;
use fs_math::dd::Dd;
use fs_ports::component::{
    CarriedStepAudit, Carrier, EntropyStepAudit, MassStepAudit, PortRole, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, audit_carried_step, audit_mass_step, audit_step,
};

use super::events::{
    EventBracket, EventDirection, EventSetConfig, EventSetError, EventSetStatus, LocatedEvent,
};
use super::{
    Checkpoint, CheckpointError, Phase, ROUNDING_UNITS, StepReport, System, SystemError,
    power_input_slack_j, power_pairing_scale,
};
use crate::solve::Condensed;
use crate::wiring::{Connection, Driver, NodeId, blocks_in_dependency_order};

/// Hard upper bound on substeps per macro window for any single rate group.
pub const MAX_MULTIRATE_SUBSTEPS: usize = 256;

/// Maximum Newton iterations allowed on a cyclic cross-rate macro interface.
pub const MAX_CROSS_RATE_ITERATIONS: u32 = 64;

/// How a slower rate group's cross-rate outputs are presented to faster groups
/// across the substeps of one macro window `[t_0, t_0 + H]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SlowInterpolation {
    /// Hold the slow group's converged macro-window output constant over every
    /// fast substep (generalizing `fs_phs::multirate::MultiratePair`). Cross-
    /// rate interface work, heat, impulse, mass, and carried quantities cancel
    /// to double-double roundoff (`interpolation_numerical_loss_j == 0.0`).
    #[default]
    ZeroOrderHold,
    /// Ramp linearly across the fast substeps from the previous accepted slow
    /// output through the macro-window mean so the window average equals the
    /// slow output identically (`u_k = y_avg + (2 * theta_k - 1) * (y_avg - y_prev)`
    /// at midpoint fraction `theta_k = (k + 0.5) / M`). On step 0 (before any
    /// prior accepted slow output exists), falls back to zero-order hold. The
    /// exact quadrature energy discrepancy `sum_k (u_k - y_avg) * y_fast_k * dt_f`
    /// is computed in double-double and booked explicitly as numerical loss.
    Linear,
}

impl SlowInterpolation {
    /// Stable canonical identifier for reports and receipts.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ZeroOrderHold => "zero_order_hold",
            Self::Linear => "linear",
        }
    }
}

/// Configuration for [`System::plan_multirate`] and [`System::step_multirate`].
#[derive(Debug, Clone, PartialEq)]
pub struct MultirateConfig {
    /// Maximum substeps allowed for any single rate group inside one macro step
    /// (`1..=MAX_MULTIRATE_SUBSTEPS`).
    pub max_substeps: usize,
    /// Round automatically inferred substep counts up to powers of two.
    pub power_of_two_ratios: bool,
    /// Interpolation law used when slower outputs feed faster groups.
    pub interpolation: SlowInterpolation,
    /// Explicit per-node substep overrides (must satisfy each node's own
    /// `limit_dt_s` and `1..=max_substeps`).
    pub node_substeps: Vec<(NodeId, usize)>,
    /// Optional explicit per-node maximum step sizes [s] (for example derived
    /// from declared modal/stiffness bounds) combined with each component's
    /// `limit_dt_s`.
    pub node_max_dt_s: Vec<(NodeId, f64)>,
    /// Maximum safeguarded Newton iterations per cyclic cross-rate macro stage
    /// (`0..=MAX_CROSS_RATE_ITERATIONS`).
    pub cross_rate_max_iterations: u32,
    /// Relative residual tolerance on cyclic cross-rate window averages.
    pub cross_rate_relative_tolerance: f64,
    /// Absolute residual tolerance on cyclic cross-rate window averages.
    pub cross_rate_absolute_tolerance: f64,
    /// When `true` (default), if a load-line loop spans nodes assigned to
    /// different rates, all nodes in that load-line block are promoted to the
    /// fastest rate among them; when `false`, splitting a load-line block is
    /// refused with [`MultirateError::SplitLoadLineBlock`].
    pub promote_load_line_blocks: bool,
}

impl Default for MultirateConfig {
    fn default() -> Self {
        Self {
            max_substeps: 64,
            power_of_two_ratios: false,
            interpolation: SlowInterpolation::ZeroOrderHold,
            node_substeps: Vec::new(),
            node_max_dt_s: Vec::new(),
            cross_rate_max_iterations: 25,
            cross_rate_relative_tolerance: 1.0e-12,
            cross_rate_absolute_tolerance: 0.0,
            promote_load_line_blocks: true,
        }
    }
}

/// Public summary of one rate group in a [`MultiratePlan`].
#[derive(Debug, Clone, PartialEq)]
pub struct RateGroup {
    /// Zero-based group index (sorted descending by `substeps`, fastest first).
    pub index: usize,
    /// Local substeps executed per macro window (`>= 1`).
    pub substeps: usize,
    /// Member components in index order.
    pub nodes: Vec<NodeId>,
    /// Member component names, parallel to `nodes`.
    pub node_names: Vec<String>,
    /// Tightest `limit_dt_s` / stiffness step bound among members at planning [s].
    pub min_limit_dt_s: f64,
    /// Ratio of `substeps` to the slowest group's `substeps` (`>= 1.0`).
    pub stiffness_ratio: f64,
}

/// Public description of one scalar output crossing between rate groups.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CrossRateChannel {
    /// Global output index inside the sealed [`System`].
    pub output_index: usize,
    /// Owning rate group index.
    pub source_group: usize,
    /// Owning component.
    pub source_node: NodeId,
    /// Owning component name.
    pub source_name: String,
    /// Downstream rate group indices that read this output.
    pub target_groups: Vec<usize>,
    /// Whether this channel is a cyclic back-edge closed by macro Newton.
    pub cyclic_feedback: bool,
}

/// Local Tarjan sub-block inside one rate group.
#[derive(Clone)]
struct LocalSubBlock {
    nodes: Vec<usize>,
    unknowns: Vec<usize>,
    known: Vec<usize>,
    /// Original `System::blocks` index when this sub-block is identical to a
    /// whole load-line block, allowing direct delegation to `System::solve_block`.
    whole_load_line_block: Option<usize>,
    /// Condensed solver pattern and output->unknown-row adjacency for looped
    /// sub-blocks without load lines.
    condensed_pattern: Vec<Vec<usize>>,
    readers: BTreeMap<usize, Vec<usize>>,
}

/// Stage in topological order across rate groups.
#[derive(Clone)]
struct GroupStage {
    /// Group indices in descending-substep (fastest-first) order.
    groups: Vec<usize>,
    /// Whether this stage has cyclic cross-rate coupling among its groups.
    cyclic: bool,
    /// Global output indices of back-edges (owned by a later/slower group in
    /// `groups` and read by an earlier/faster group in `groups`).
    back_outputs: Vec<usize>,
}

/// Precomputed multirate schedule for one sealed [`System`].
#[derive(Clone)]
pub struct MultiratePlan {
    identity: Rc<()>,
    config: MultirateConfig,
    groups: Vec<RateGroup>,
    node_group: Vec<usize>,
    group_sub_blocks: Vec<Vec<LocalSubBlock>>,
    group_boundaries: Vec<Vec<usize>>,
    stages: Vec<GroupStage>,
    cross_channels: Vec<CrossRateChannel>,
    /// Per global input index, whether its driver reads any output owned by a
    /// strictly slower group (`substeps_src < substeps_dst`).
    input_has_slower_remote: Vec<bool>,
    /// Per global output index, whether it is ever read by a temperature input.
    output_feeds_temperature: Vec<bool>,
}

impl fmt::Debug for MultiratePlan {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MultiratePlan")
            .field("groups", &self.groups)
            .field("cross_channels", &self.cross_channels)
            .field("interpolation", &self.config.interpolation)
            .field("stages", &self.stages.len())
            .finish_non_exhaustive()
    }
}

impl MultiratePlan {
    /// Whether every rate group runs at `substeps == 1` (single-rate).
    #[must_use]
    pub fn is_single_rate(&self) -> bool {
        self.groups.iter().all(|g| g.substeps == 1)
    }

    /// Rate groups in descending substep order (fastest first).
    #[must_use]
    pub fn groups(&self) -> &[RateGroup] {
        &self.groups
    }

    /// Cross-rate scalar output channels between different rate groups.
    #[must_use]
    pub fn cross_rate_channels(&self) -> &[CrossRateChannel] {
        &self.cross_channels
    }

    /// Declared slow-to-fast interpolation law.
    #[must_use]
    pub const fn interpolation(&self) -> SlowInterpolation {
        self.config.interpolation
    }

    /// Maximum substeps across all rate groups.
    #[must_use]
    pub fn max_substeps(&self) -> usize {
        self.groups.iter().map(|g| g.substeps).max().unwrap_or(1)
    }

    /// Number of cyclic cross-rate stages requiring macro interface iteration.
    #[must_use]
    pub fn cyclic_stage_count(&self) -> usize {
        self.stages.iter().filter(|s| s.cyclic).count()
    }

    /// Rate group index assigned to `node`, if `node` belongs to this system.
    #[must_use]
    pub fn node_group(&self, node: NodeId) -> Option<usize> {
        self.node_group.get(node.0).copied()
    }
}

/// Summary of one accepted multirate macro step.
#[derive(Debug, Clone, PartialEq)]
pub struct MultirateStepReport {
    /// Accepted macro step duration `H` [s].
    pub dt_s: f64,
    /// Accepted absolute simulation time after publication [s].
    pub time_s: f64,
    /// Substeps executed by each rate group (matching [`MultiratePlan::groups`]).
    pub group_substeps: Vec<usize>,
    /// Total cyclic cross-rate macro Newton iterations spent across stages.
    pub cross_rate_iterations: u32,
    /// Total internal sub-block Newton iterations spent across all substeps.
    pub inner_iterations: u32,
    /// Macro shortenings performed for component/boundary events or requests.
    pub shortenings: u32,
    /// Exact energy numerical loss booked from slow-to-fast interpolation [J].
    /// Identically `0.0` under [`SlowInterpolation::ZeroOrderHold`].
    pub interpolation_numerical_loss_j: f64,
    /// Exact entropy discrepancy from slow-to-fast heat-port temperature
    /// interpolation [J/K]. Identically `0.0` under `ZeroOrderHold`.
    pub interpolation_entropy_defect_j_per_k: f64,
    /// Largest scaled cross-rate residual on convergence (`0.0` if acyclic).
    pub max_cross_rate_residual: f64,
}

/// Accepted result of [`System::step_multirate_to_events`].
#[derive(Debug, Clone, PartialEq)]
pub struct MultirateEventReport {
    /// Whether any guard crossed zero during the macro window.
    pub status: EventSetStatus,
    /// The single retained synchronized multirate macro step, or `None` on
    /// entry zeros.
    pub step: Option<MultirateStepReport>,
    /// Accepted synchronized simulation time shared by every rate group [s].
    pub time_s: f64,
    /// Every guard's value at the accepted state, in declaration order.
    pub values: Vec<f64>,
    /// Completed trial macro steps, shared by all guards.
    pub trial_steps: u64,
}

/// Why a multirate plan or step was refused. Every step refusal restores the
/// exact entry state, clock, and ledger.
#[derive(Debug)]
pub enum MultirateError {
    /// An invalid multirate configuration field or node index.
    InvalidConfig {
        /// Offending field name or detail.
        field: &'static str,
    },
    /// A component's `limit_dt_s` (or `node_max_dt_s`) requires more substeps
    /// than `config.max_substeps`.
    RateRatioExceeded {
        /// Component name.
        node: String,
        /// Required substeps for the macro window.
        required: usize,
        /// Configured ceiling.
        max_substeps: usize,
    },
    /// An explicit `node_substeps` override is coarser than the component's own
    /// `limit_dt_s` allows.
    OverrideViolatesLimit {
        /// Component name.
        node: String,
        /// Requested override substeps.
        requested: usize,
        /// Minimum substeps required by `limit_dt_s`.
        required: usize,
    },
    /// A load-line block was split across distinct rate groups with
    /// `promote_load_line_blocks == false`.
    SplitLoadLineBlock {
        /// Component names in the load-line block.
        nodes: Vec<String>,
    },
    /// The supplied [`MultiratePlan`] belongs to a different [`System`].
    ForeignPlan,
    /// Slow-to-fast linear interpolation produced a nonfinite or non-positive
    /// temperature value.
    InvalidInterpolation {
        /// Output owner name.
        node: String,
        /// Interpolated value that was refused.
        value: f64,
    },
    /// A cyclic cross-rate macro interface failed to converge within budget.
    CrossRateNotConverged {
        /// Iterations spent.
        iterations: u32,
        /// Largest relative residual remaining.
        residual: f64,
    },
    /// A cyclic cross-rate macro Jacobian was singular or nonfinite.
    CrossRateSingular,
    /// Underlying kernel step/audit/negotiation refusal.
    System(Box<SystemError>),
    /// Event-set configuration or bisection refusal.
    Event(Box<EventSetError>),
    /// Checkpoint capture or restore refusal.
    Checkpoint(Box<CheckpointError>),
}

impl fmt::Display for MultirateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field } => {
                write!(f, "invalid multirate configuration: {field}")
            }
            Self::RateRatioExceeded {
                node,
                required,
                max_substeps,
            } => write!(
                f,
                "multirate rate ratio exceeded for `{node}`: requires {required} substeps (max {max_substeps})"
            ),
            Self::OverrideViolatesLimit {
                node,
                requested,
                required,
            } => write!(
                f,
                "multirate override for `{node}` ({requested} substeps) is coarser than its limit_dt_s requires ({required} substeps)"
            ),
            Self::SplitLoadLineBlock { nodes } => write!(
                f,
                "multirate plan splits load-line block {nodes:?} across different rate groups"
            ),
            Self::ForeignPlan => write!(f, "multirate plan belongs to a different System"),
            Self::InvalidInterpolation { node, value } => write!(
                f,
                "multirate linear interpolation for `{node}` produced invalid value {value}"
            ),
            Self::CrossRateNotConverged {
                iterations,
                residual,
            } => write!(
                f,
                "cross-rate macro interface did not converge in {iterations} iterations (residual {residual:.3e})"
            ),
            Self::CrossRateSingular => write!(f, "cross-rate macro interface Jacobian is singular"),
            Self::System(err) => fmt::Display::fmt(err, f),
            Self::Event(err) => fmt::Display::fmt(err, f),
            Self::Checkpoint(err) => fmt::Display::fmt(err, f),
        }
    }
}

impl std::error::Error for MultirateError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::System(err) => Some(err.as_ref()),
            Self::Event(err) => Some(err.as_ref()),
            Self::Checkpoint(err) => Some(err.as_ref()),
            _ => None,
        }
    }
}

impl From<SystemError> for MultirateError {
    fn from(err: SystemError) -> Self {
        Self::System(Box::new(err))
    }
}

impl From<EventSetError> for MultirateError {
    fn from(err: EventSetError) -> Self {
        Self::Event(Box::new(err))
    }
}

impl From<CheckpointError> for MultirateError {
    fn from(err: CheckpointError) -> Self {
        Self::Checkpoint(Box::new(err))
    }
}

/// Double-double accumulators for one node across all its substeps in a macro
/// window.
#[derive(Clone, Default)]
struct NodeMacroAccum {
    port_energy_in_j: Dd,
    stored_change_j: Dd,
    dissipated_j: Dd,
    numerical_loss_j: Dd,
    entropy_generated_j_per_k: Dd,
    balance_defect_j: Dd,
    absolute_defect_j: Dd,
    interp_energy_j: Dd,
    has_mass: bool,
    mass_port_in_kg: Dd,
    mass_stored_change_kg: Dd,
    mass_balance_defect_kg: Dd,
    mass_abs_defect_kg: Dd,
    mass_audited_steps: u64,
    has_entropy: bool,
    entropy_port_in: Dd,
    entropy_gen_audit: Dd,
    entropy_stored_change: Dd,
    entropy_balance_defect: Dd,
    entropy_abs_defect: Dd,
    entropy_audited_steps: u64,
    carried: Vec<Option<CarriedMacroAccum>>,
}

#[derive(Clone, Default)]
struct CarriedMacroAccum {
    port_in: Dd,
    stored_change: Dd,
    balance_defect: Dd,
    abs_defect: Dd,
    audited_steps: u64,
}

/// Double-double accumulators collected across one complete macro-window trial.
struct MacroAccum {
    nodes: Vec<NodeMacroAccum>,
    port_energy: Vec<Dd>,
    port_entropy: Vec<Dd>,
    port_mass: Vec<Dd>,
    port_carried: Vec<Vec<Dd>>,
    port_magnitude: Vec<Dd>,
    port_interp_energy: Vec<Dd>,
    port_interp_entropy: Vec<Dd>,
    port_solve_slack_j: Vec<f64>,
    port_mass_slack_kg: Vec<f64>,
    port_entropy_slack: Vec<f64>,
    port_carried_slack: Vec<Vec<f64>>,
    window_outputs: Vec<f64>,
    window_inputs: Vec<f64>,
    window_allowance: Vec<f64>,
    inner_iterations: u32,
}

impl MacroAccum {
    fn new(system: &System) -> Self {
        let n_nodes = system.nodes.len();
        let n_ports = system.port_energy.len();
        let mut nodes = Vec::with_capacity(n_nodes);
        for carried in &system.carried_audits {
            nodes.push(NodeMacroAccum {
                carried: vec![None; carried.len()],
                ..NodeMacroAccum::default()
            });
        }
        let port_carried = system
            .port_carried
            .iter()
            .map(|c| vec![Dd::ZERO; c.len()])
            .collect();
        let port_carried_slack = system
            .port_carried
            .iter()
            .map(|c| vec![0.0_f64; c.len()])
            .collect();
        Self {
            nodes,
            port_energy: vec![Dd::ZERO; n_ports],
            port_entropy: vec![Dd::ZERO; n_ports],
            port_mass: vec![Dd::ZERO; n_ports],
            port_carried,
            port_magnitude: vec![Dd::ZERO; n_ports],
            port_interp_energy: vec![Dd::ZERO; n_ports],
            port_interp_entropy: vec![Dd::ZERO; n_ports],
            port_solve_slack_j: vec![0.0; n_ports],
            port_mass_slack_kg: vec![0.0; n_ports],
            port_entropy_slack: vec![0.0; n_ports],
            port_carried_slack,
            window_outputs: system.accepted_outputs.clone(),
            window_inputs: system.accepted_inputs.clone(),
            window_allowance: vec![0.0; system.inputs.len()],
            inner_iterations: 0,
        }
    }

    fn reset_group(&mut self, system: &System, nodes: &[usize]) {
        for &n in nodes {
            self.nodes[n] = NodeMacroAccum {
                carried: vec![None; system.carried_audits[n].len()],
                ..NodeMacroAccum::default()
            };
            for global in system.port_start[n]..system.port_start[n + 1] {
                self.port_energy[global] = Dd::ZERO;
                self.port_entropy[global] = Dd::ZERO;
                self.port_mass[global] = Dd::ZERO;
                self.port_carried[global].fill(Dd::ZERO);
                self.port_magnitude[global] = Dd::ZERO;
                self.port_interp_energy[global] = Dd::ZERO;
                self.port_interp_entropy[global] = Dd::ZERO;
                self.port_solve_slack_j[global] = 0.0;
                self.port_mass_slack_kg[global] = 0.0;
                self.port_entropy_slack[global] = 0.0;
                self.port_carried_slack[global].fill(0.0);
            }
        }
    }
}

enum GroupOutcome {
    Done,
    Shorten { macro_dt_s: f64 },
}

enum SubBlockOutcome {
    Done { iterations: u32 },
    Shorten { dt_s: f64 },
}

fn validate_config(config: &MultirateConfig, node_count: usize) -> Result<(), MultirateError> {
    if !(1..=MAX_MULTIRATE_SUBSTEPS).contains(&config.max_substeps) {
        return Err(MultirateError::InvalidConfig {
            field: "max_substeps",
        });
    }
    if config.cross_rate_max_iterations > MAX_CROSS_RATE_ITERATIONS {
        return Err(MultirateError::InvalidConfig {
            field: "cross_rate_max_iterations",
        });
    }
    if !(config.cross_rate_relative_tolerance.is_finite()
        && config.cross_rate_relative_tolerance > 0.0)
    {
        return Err(MultirateError::InvalidConfig {
            field: "cross_rate_relative_tolerance",
        });
    }
    if !(config.cross_rate_absolute_tolerance.is_finite()
        && config.cross_rate_absolute_tolerance >= 0.0)
    {
        return Err(MultirateError::InvalidConfig {
            field: "cross_rate_absolute_tolerance",
        });
    }
    for &(node, substeps) in &config.node_substeps {
        if node.0 >= node_count || !(1..=config.max_substeps).contains(&substeps) {
            return Err(MultirateError::InvalidConfig {
                field: "node_substeps",
            });
        }
    }
    for &(node, dt) in &config.node_max_dt_s {
        if node.0 >= node_count || !(dt.is_finite() && dt > 0.0) {
            return Err(MultirateError::InvalidConfig {
                field: "node_max_dt_s",
            });
        }
    }
    Ok(())
}

impl System {
    /// Assign components to rate groups for a macro step of `macro_dt_s`,
    /// combining each component's `limit_dt_s` with any explicit stiffness /
    /// substep overrides in `config`.
    ///
    /// # Errors
    /// Invalid duration or config, component negotiation refusal, exceeded rate
    /// ratio ceiling, override coarser than a component's own `limit_dt_s`, or
    /// a split load-line block when `promote_load_line_blocks` is `false`.
    #[allow(clippy::too_many_lines)]
    pub fn plan_multirate(
        &self,
        macro_dt_s: f64,
        config: MultirateConfig,
    ) -> Result<MultiratePlan, MultirateError> {
        if !(macro_dt_s.is_finite() && macro_dt_s > 0.0) {
            return Err(SystemError::Step { dt_s: macro_dt_s }.into());
        }
        let n_nodes = self.nodes.len();
        validate_config(&config, n_nodes)?;

        let mut max_dt_overrides: BTreeMap<usize, f64> = BTreeMap::new();
        for &(node, dt) in &config.node_max_dt_s {
            max_dt_overrides
                .entry(node.0)
                .and_modify(|cur| *cur = cur.min(dt))
                .or_insert(dt);
        }
        let mut substep_overrides: BTreeMap<usize, usize> = BTreeMap::new();
        for &(node, substeps) in &config.node_substeps {
            substep_overrides.insert(node.0, substeps);
        }

        let mut node_substeps = vec![1_usize; n_nodes];
        let mut node_limits = vec![macro_dt_s; n_nodes];
        for node in 0..n_nodes {
            let limited = self.nodes[node]
                .limit_dt_s(self.accepted[node].as_ref(), macro_dt_s)
                .map_err(|error| SystemError::Component {
                    node: self.names[node].clone(),
                    phase: Phase::Negotiate,
                    error,
                })?;
            if !(limited > 0.0 && limited <= macro_dt_s) {
                return Err(SystemError::Negotiation {
                    node: self.names[node].clone(),
                    proposed_s: macro_dt_s,
                    returned_s: limited,
                }
                .into());
            }
            let eff_dt = max_dt_overrides
                .get(&node)
                .copied()
                .map_or(limited, |cap| limited.min(cap));
            node_limits[node] = eff_dt;

            let ratio = macro_dt_s / eff_dt;
            if !ratio.is_finite() {
                return Err(MultirateError::RateRatioExceeded {
                    node: self.names[node].clone(),
                    required: usize::MAX,
                    max_substeps: config.max_substeps,
                });
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let natural = if ratio <= 1.0 + 1.0e-12 {
                1_usize
            } else if ratio.ceil() > config.max_substeps as f64 {
                return Err(MultirateError::RateRatioExceeded {
                    node: self.names[node].clone(),
                    required: ratio.ceil() as usize,
                    max_substeps: config.max_substeps,
                });
            } else {
                ratio.ceil() as usize
            };

            let mut chosen = if let Some(&requested) = substep_overrides.get(&node) {
                if requested < natural {
                    return Err(MultirateError::OverrideViolatesLimit {
                        node: self.names[node].clone(),
                        requested,
                        required: natural,
                    });
                }
                requested
            } else {
                natural
            };
            if config.power_of_two_ratios {
                chosen = chosen.next_power_of_two();
                if chosen > config.max_substeps {
                    return Err(MultirateError::RateRatioExceeded {
                        node: self.names[node].clone(),
                        required: chosen,
                        max_substeps: config.max_substeps,
                    });
                }
            }
            node_substeps[node] = chosen;
        }

        // Keep every load-line loop intact in one rate group.
        for (b_idx, block) in self.blocks.iter().enumerate() {
            if self.lines[b_idx].is_empty() || block.nodes.len() <= 1 {
                continue;
            }
            let min_s = block
                .nodes
                .iter()
                .map(|&n| node_substeps[n])
                .min()
                .unwrap_or(1);
            let max_s = block
                .nodes
                .iter()
                .map(|&n| node_substeps[n])
                .max()
                .unwrap_or(1);
            if min_s != max_s {
                if config.promote_load_line_blocks {
                    for &n in &block.nodes {
                        node_substeps[n] = max_s;
                    }
                } else {
                    return Err(MultirateError::SplitLoadLineBlock {
                        nodes: self.names_of(&block.nodes),
                    });
                }
            }
        }

        // Collect distinct rates, fastest first.
        let mut distinct_rates: Vec<usize> = node_substeps.clone();
        distinct_rates.sort_unstable();
        distinct_rates.dedup();
        distinct_rates.reverse();

        let slowest_substeps = *distinct_rates.last().unwrap_or(&1) as f64;
        let mut node_group = vec![0_usize; n_nodes];
        let mut groups = Vec::with_capacity(distinct_rates.len());
        for (g_idx, &substeps) in distinct_rates.iter().enumerate() {
            let mut g_nodes = Vec::new();
            let mut min_limit = f64::INFINITY;
            for (n, &s) in node_substeps.iter().enumerate() {
                if s == substeps {
                    node_group[n] = g_idx;
                    g_nodes.push(NodeId(n));
                    min_limit = min_limit.min(node_limits[n]);
                }
            }
            let node_names = g_nodes.iter().map(|id| self.names[id.0].clone()).collect();
            groups.push(RateGroup {
                index: g_idx,
                substeps,
                nodes: g_nodes,
                node_names,
                min_limit_dt_s: min_limit,
                stiffness_ratio: (substeps as f64) / slowest_substeps,
            });
        }

        let output_owner = self.multirate_output_owners();
        let output_group: Vec<usize> = output_owner.iter().map(|&n| node_group[n]).collect();

        let group_sub_blocks = self.build_group_sub_blocks(&groups, &node_group, &output_owner);
        let group_boundaries = self.build_group_boundaries(groups.len(), &node_group);

        // Identify cross-rate edges and build topological stages across groups.
        let (stages, cross_channels, input_has_slower_remote, output_feeds_temperature) =
            self.build_cross_rate_topology(&groups, &node_group, &output_owner, &output_group);

        Ok(MultiratePlan {
            identity: Rc::clone(&self.identity),
            config,
            groups,
            node_group,
            group_sub_blocks,
            group_boundaries,
            stages,
            cross_channels,
            input_has_slower_remote,
            output_feeds_temperature,
        })
    }

    fn multirate_output_owners(&self) -> Vec<usize> {
        let mut owners = Vec::with_capacity(self.outputs.len());
        for node in 0..self.nodes.len() {
            owners.extend(core::iter::repeat_n(
                node,
                self.output_start[node + 1] - self.output_start[node],
            ));
        }
        owners
    }

    fn build_group_boundaries(&self, group_count: usize, node_group: &[usize]) -> Vec<Vec<usize>> {
        let mut out = vec![Vec::new(); group_count];
        for (b_idx, (port, _)) in self.boundaries.iter().enumerate() {
            out[node_group[port.node]].push(b_idx);
        }
        out
    }

    fn build_group_sub_blocks(
        &self,
        groups: &[RateGroup],
        node_group: &[usize],
        output_owner: &[usize],
    ) -> Vec<Vec<LocalSubBlock>> {
        let mut result = Vec::with_capacity(groups.len());
        let mut src_buf = Vec::new();
        for (g_idx, group) in groups.iter().enumerate() {
            let local_nodes: Vec<usize> = group.nodes.iter().map(|id| id.0).collect();
            let mut local_pos = BTreeMap::new();
            for (pos, &node) in local_nodes.iter().enumerate() {
                local_pos.insert(node, pos);
            }
            let mut successors: Vec<Vec<usize>> = vec![Vec::new(); local_nodes.len()];
            for (dst_pos, &dst_node) in local_nodes.iter().enumerate() {
                for in_idx in self.node_inputs(dst_node) {
                    src_buf.clear();
                    self.drivers[in_idx].sources(&mut src_buf);
                    for &out_idx in &src_buf {
                        let src_node = output_owner[out_idx];
                        if node_group[src_node] == g_idx
                            && let Some(&src_pos) = local_pos.get(&src_node)
                            && !successors[src_pos].contains(&dst_pos)
                        {
                            successors[src_pos].push(dst_pos);
                        }
                    }
                }
            }
            let sccs = blocks_in_dependency_order(&successors);
            let mut sub_blocks = Vec::with_capacity(sccs.len());
            for scc in sccs {
                let mut nodes: Vec<usize> = scc.into_iter().map(|p| local_nodes[p]).collect();
                nodes.sort_unstable();
                sub_blocks.push(self.build_one_sub_block(nodes, output_owner, &mut src_buf));
            }
            result.push(sub_blocks);
        }
        result
    }

    fn build_one_sub_block(
        &self,
        nodes: Vec<usize>,
        output_owner: &[usize],
        src_buf: &mut Vec<usize>,
    ) -> LocalSubBlock {
        let mut unknowns = Vec::new();
        let mut known = Vec::new();
        for &node in &nodes {
            for in_idx in self.node_inputs(node) {
                src_buf.clear();
                self.drivers[in_idx].sources(src_buf);
                let depends_in_block = src_buf
                    .iter()
                    .any(|&out_idx| nodes.binary_search(&output_owner[out_idx]).is_ok());
                if depends_in_block {
                    unknowns.push(in_idx);
                } else {
                    known.push(in_idx);
                }
            }
        }
        let whole_load_line_block = self
            .blocks
            .iter()
            .enumerate()
            .find(|(b_idx, b)| !self.lines[*b_idx].is_empty() && b.nodes == nodes)
            .map(|(b_idx, _)| b_idx);

        let mut readers: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut pattern = vec![Vec::new(); unknowns.len()];
        if !unknowns.is_empty() && whole_load_line_block.is_none() {
            for (row, &in_idx) in unknowns.iter().enumerate() {
                src_buf.clear();
                self.drivers[in_idx].sources(src_buf);
                for &out_idx in src_buf.iter() {
                    if nodes.binary_search(&output_owner[out_idx]).is_ok() {
                        readers.entry(out_idx).or_default().push(row);
                    }
                }
            }
            for (col, &in_idx) in unknowns.iter().enumerate() {
                let owner = self.input_owner[in_idx];
                for out_idx in self.node_outputs(owner) {
                    if let Some(rows) = readers.get(&out_idx) {
                        for &row in rows {
                            if !pattern[col].contains(&row) {
                                pattern[col].push(row);
                            }
                        }
                    }
                }
                pattern[col].sort_unstable();
            }
        }

        LocalSubBlock {
            nodes,
            unknowns,
            known,
            whole_load_line_block,
            condensed_pattern: pattern,
            readers,
        }
    }

    fn build_cross_rate_topology(
        &self,
        groups: &[RateGroup],
        node_group: &[usize],
        output_owner: &[usize],
        output_group: &[usize],
    ) -> (Vec<GroupStage>, Vec<CrossRateChannel>, Vec<bool>, Vec<bool>) {
        let g_count = groups.len();
        let mut g_succ: Vec<Vec<usize>> = vec![Vec::new(); g_count];
        let mut cross_targets: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut input_has_slower_remote = vec![false; self.inputs.len()];
        let mut output_feeds_temperature = vec![false; self.outputs.len()];
        let mut src_buf = Vec::new();

        for (in_idx, driver) in self.drivers.iter().enumerate() {
            let dst_node = self.input_owner[in_idx];
            let dst_group = node_group[dst_node];
            let is_temp = self.input_is_temperature[in_idx];
            src_buf.clear();
            driver.sources(&mut src_buf);
            for &out_idx in &src_buf {
                if is_temp {
                    output_feeds_temperature[out_idx] = true;
                }
                let src_group = output_group[out_idx];
                if src_group != dst_group {
                    if !g_succ[src_group].contains(&dst_group) {
                        g_succ[src_group].push(dst_group);
                    }
                    let targets = cross_targets.entry(out_idx).or_default();
                    if !targets.contains(&dst_group) {
                        targets.push(dst_group);
                    }
                    if groups[src_group].substeps < groups[dst_group].substeps {
                        input_has_slower_remote[in_idx] = true;
                    }
                }
            }
        }

        let sccs = blocks_in_dependency_order(&g_succ);
        let mut stages = Vec::with_capacity(sccs.len());
        let mut back_set = Vec::new();
        for mut scc in sccs {
            // Sort fastest first (lowest group index has highest `substeps`).
            scc.sort_unstable();
            let cyclic = scc.len() > 1;
            let mut back_outputs = Vec::new();
            if cyclic {
                for (&out_idx, targets) in &cross_targets {
                    let src_g = output_group[out_idx];
                    if let Some(src_pos) = scc.iter().position(|&g| g == src_g) {
                        let feeds_earlier = targets.iter().any(|&tg| {
                            scc.iter()
                                .position(|&g| g == tg)
                                .is_some_and(|t_pos| t_pos < src_pos)
                        });
                        if feeds_earlier {
                            back_outputs.push(out_idx);
                            back_set.push(out_idx);
                        }
                    }
                }
            }
            stages.push(GroupStage {
                groups: scc,
                cyclic,
                back_outputs,
            });
        }

        let cross_channels = cross_targets
            .into_iter()
            .map(|(output_index, mut target_groups)| {
                target_groups.sort_unstable();
                let source_node = output_owner[output_index];
                CrossRateChannel {
                    output_index,
                    source_group: output_group[output_index],
                    source_node: NodeId(source_node),
                    source_name: self.names[source_node].clone(),
                    target_groups,
                    cyclic_feedback: back_set.contains(&output_index),
                }
            })
            .collect();

        (
            stages,
            cross_channels,
            input_has_slower_remote,
            output_feeds_temperature,
        )
    }

    /// Advance the system by at most `macro_dt_s` using the multirate schedule
    /// `plan`.
    ///
    /// When `plan.is_single_rate()` is `true`, this delegates directly to
    /// [`Self::step`], reproducing single-rate evolution bit-for-bit.
    ///
    /// Every refusal or cancellation restores the exact accepted state, clock,
    /// and ledger from entry (zero partial commit).
    ///
    /// # Errors
    /// Foreign plan, invalid duration, component/interface/second-law refusal,
    /// cross-rate non-convergence, or cancellation.
    pub fn step_multirate(
        &mut self,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        cx: &Cx<'_>,
    ) -> Result<MultirateStepReport, MultirateError> {
        if !Rc::ptr_eq(&self.identity, &plan.identity) {
            return Err(MultirateError::ForeignPlan);
        }
        if !(macro_dt_s.is_finite() && macro_dt_s > 0.0) {
            return Err(SystemError::Step { dt_s: macro_dt_s }.into());
        }
        if plan.is_single_rate() {
            let report = self.step(macro_dt_s, cx)?;
            return Ok(MultirateStepReport {
                dt_s: report.dt_s,
                time_s: report.time_s,
                group_substeps: vec![1; plan.groups.len()],
                cross_rate_iterations: 0,
                inner_iterations: report.iterations,
                shortenings: report.shortenings,
                interpolation_numerical_loss_j: 0.0,
                interpolation_entropy_defect_j_per_k: 0.0,
                max_cross_rate_residual: 0.0,
            });
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;

        let mut h = self.negotiate_macro_window(macro_dt_s, plan)?;
        let entry = self.checkpoint()?;
        let entry_allowance = self.accepted_allowance.clone();
        let entry_receipts = self.accepted_receipts.clone();

        let result = self.run_multirate_trials(&mut h, plan, &entry, cx);
        if result.is_err() {
            self.restore(&entry)?;
            self.accepted_allowance.copy_from_slice(&entry_allowance);
            self.accepted_receipts.copy_from_slice(&entry_receipts);
        }
        result
    }

    fn negotiate_macro_window(
        &self,
        macro_dt_s: f64,
        plan: &MultiratePlan,
    ) -> Result<f64, MultirateError> {
        let mut h = macro_dt_s;
        for (g_idx, group) in plan.groups.iter().enumerate() {
            let dt_g = h / (group.substeps as f64);
            for &node_id in &group.nodes {
                let node = node_id.0;
                let limited = self.nodes[node]
                    .limit_dt_s(self.accepted[node].as_ref(), dt_g)
                    .map_err(|error| SystemError::Component {
                        node: self.names[node].clone(),
                        phase: Phase::Negotiate,
                        error,
                    })?;
                if !(limited > 0.0 && limited <= dt_g) {
                    return Err(SystemError::Negotiation {
                        node: self.names[node].clone(),
                        proposed_s: dt_g,
                        returned_s: limited,
                    }
                    .into());
                }
                if limited < dt_g {
                    h = h.min(limited * (group.substeps as f64));
                }
                if let Some(event_s) = self.nodes[node].next_event_s(self.accepted[node].as_ref()) {
                    if !(event_s > 0.0) {
                        return Err(SystemError::Negotiation {
                            node: self.names[node].clone(),
                            proposed_s: h,
                            returned_s: event_s,
                        }
                        .into());
                    }
                    h = h.min(event_s);
                }
            }
            let _ = g_idx;
        }
        // Clip to the earliest scheduled boundary event across ALL groups so
        // every rate group synchronizes on the exact event clock.
        let next_boundary = self
            .boundary_events
            .iter()
            .enumerate()
            .filter_map(|(boundary, events)| {
                let next = events.partition_point(|&event| event <= self.time.value());
                events.get(next).map(|&time_s| (boundary, time_s))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1));
        if let Some((boundary, time_s)) = next_boundary {
            let event_s = self.time.interval_to(time_s);
            if !(event_s > 0.0) {
                return Err(SystemError::Negotiation {
                    node: self.ledger.boundaries[boundary].port.clone(),
                    proposed_s: h,
                    returned_s: event_s,
                }
                .into());
            }
            h = self.time.limit_before_event(h, event_s).min(event_s);
        }
        Ok(h)
    }

    fn run_multirate_trials(
        &mut self,
        h: &mut f64,
        plan: &MultiratePlan,
        entry: &Checkpoint,
        cx: &Cx<'_>,
    ) -> Result<MultirateStepReport, MultirateError> {
        let mut shortenings = 0_u32;
        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            let mut accum = MacroAccum::new(self);
            let mut cross_iters = 0_u32;
            let mut max_cross_res = 0.0_f64;
            let mut shortened_to = None;

            for stage in &plan.stages {
                if stage.cyclic {
                    match self.solve_cyclic_stage(stage, *h, plan, entry, &mut accum, cx)? {
                        CyclicStageOutcome::Done {
                            iterations,
                            max_residual,
                        } => {
                            cross_iters += iterations;
                            max_cross_res = max_cross_res.max(max_residual);
                        }
                        CyclicStageOutcome::Shorten { macro_dt_s } => {
                            shortened_to = Some(macro_dt_s);
                            break;
                        }
                    }
                } else {
                    let g = stage.groups[0];
                    match self.eval_rate_group(g, *h, plan, entry, &mut accum, cx)? {
                        GroupOutcome::Done => {}
                        GroupOutcome::Shorten { macro_dt_s } => {
                            shortened_to = Some(macro_dt_s);
                            break;
                        }
                    }
                }
            }

            if let Some(shorter_h) = shortened_to {
                shortenings += 1;
                if shortenings > self.config.max_shortenings || !(shorter_h > 0.0 && shorter_h < *h)
                {
                    return Err(SystemError::StepUnderflow { dt_s: shorter_h }.into());
                }
                self.restore(entry)?;
                *h = shorter_h;
                continue;
            }

            return self.audit_and_publish_macro(
                *h,
                plan,
                &accum,
                cross_iters,
                shortenings,
                max_cross_res,
                cx,
            );
        }
    }

    fn restore_group_nodes(
        &mut self,
        nodes: &[usize],
        entry: &Checkpoint,
    ) -> Result<(), MultirateError> {
        for &n in nodes {
            self.accepted[n] = self.nodes[n]
                .clone_state(entry.accepted[n].as_ref())
                .map_err(|source| CheckpointError::State {
                    node: self.names[n].clone(),
                    source,
                })?;
            let in_r = self.node_inputs(n);
            self.inputs[in_r.clone()].copy_from_slice(&entry.inputs[in_r]);
            let out_r = self.node_outputs(n);
            self.outputs[out_r.clone()].copy_from_slice(&entry.outputs[out_r]);
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn eval_rate_group(
        &mut self,
        group_idx: usize,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        entry: &Checkpoint,
        accum: &mut MacroAccum,
        cx: &Cx<'_>,
    ) -> Result<GroupOutcome, MultirateError> {
        let group = &plan.groups[group_idx];
        let g_nodes: Vec<usize> = group.nodes.iter().map(|id| id.0).collect();
        self.restore_group_nodes(&g_nodes, entry)?;
        accum.reset_group(self, &g_nodes);

        let substeps = group.substeps;
        let dt_g = macro_dt_s / (substeps as f64);
        let mut t_sub = entry.time;

        let mut output_integrals: BTreeMap<usize, Dd> = BTreeMap::new();
        let mut input_integrals: BTreeMap<usize, Dd> = BTreeMap::new();
        for &n in &g_nodes {
            for o in self.node_outputs(n) {
                output_integrals.insert(o, Dd::ZERO);
            }
            for i in self.node_inputs(n) {
                input_integrals.insert(i, Dd::ZERO);
            }
        }

        let mut delta_remote = vec![0.0_f64; self.outputs.len()];
        let use_linear = plan.config.interpolation == SlowInterpolation::Linear
            && entry.steps > 0
            && substeps > 1;

        for k in 0..substeps {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            // 1. Evaluate boundaries owned by this rate group at its local clock.
            for &b_idx in &plan.group_boundaries[group_idx] {
                let t_now = t_sub.value();
                (self.boundaries[b_idx].1)(t_now, dt_g, &mut self.boundary_values[b_idx]);
                for (v, &offset) in self.boundary_values[b_idx]
                    .iter_mut()
                    .zip(&self.boundary_bias[b_idx])
                {
                    if offset != 0.0 {
                        *v += offset;
                    }
                }
            }

            // 2. Populate remote outputs in `self.outputs` for substep `k`.
            delta_remote.fill(0.0);
            let theta_k = ((k as f64) + 0.5) / (substeps as f64);
            let ramp = 2.0 * theta_k - 1.0;
            for ch in &plan.cross_channels {
                if !ch.target_groups.contains(&group_idx) {
                    continue;
                }
                let o = ch.output_index;
                let y_avg = accum.window_outputs[o];
                let src_substeps = plan.groups[ch.source_group].substeps;
                if use_linear && src_substeps < substeps {
                    let y_prev = entry.outputs[o];
                    let delta = ramp * (y_avg - y_prev);
                    let val = y_avg + delta;
                    if !val.is_finite() || (plan.output_feeds_temperature[o] && !(val > 0.0)) {
                        return Err(MultirateError::InvalidInterpolation {
                            node: ch.source_name.clone(),
                            value: val,
                        });
                    }
                    self.outputs[o] = val;
                    delta_remote[o] = delta;
                } else {
                    self.outputs[o] = y_avg;
                }
            }

            // Clear allowances for this group's inputs before local sub-blocks.
            for &n in &g_nodes {
                let in_r = self.node_inputs(n);
                self.input_allowance[in_r].fill(0.0);
            }

            // 3. Step every local sub-block in topological order.
            for sub_block in &plan.group_sub_blocks[group_idx] {
                for &in_idx in &sub_block.known {
                    self.inputs[in_idx] =
                        self.drivers[in_idx].value(&self.outputs, &self.boundary_values);
                }
                if sub_block.unknowns.is_empty() {
                    for &node in &sub_block.nodes {
                        if let Some(shorter) = self.propose(node, dt_g, cx)? {
                            let elapsed = (k as f64) * dt_g + shorter;
                            return Ok(GroupOutcome::Shorten {
                                macro_dt_s: elapsed,
                            });
                        }
                    }
                } else if let Some(orig_block) = sub_block.whole_load_line_block {
                    match self.solve_block(orig_block, dt_g, cx)? {
                        Ok(iters) => accum.inner_iterations += iters,
                        Err(shorter) => {
                            let elapsed = (k as f64) * dt_g + shorter;
                            return Ok(GroupOutcome::Shorten {
                                macro_dt_s: elapsed,
                            });
                        }
                    }
                } else {
                    match self.solve_local_sub_block(sub_block, dt_g, cx)? {
                        SubBlockOutcome::Done { iterations } => {
                            accum.inner_iterations += iterations;
                        }
                        SubBlockOutcome::Shorten { dt_s: shorter } => {
                            let elapsed = (k as f64) * dt_g + shorter;
                            return Ok(GroupOutcome::Shorten {
                                macro_dt_s: elapsed,
                            });
                        }
                    }
                }
            }

            // 4. Audit & admit each node in this group for substep `k`, and
            //    accumulate exact double-double totals.
            self.audit_admit_and_accumulate_substep(
                &g_nodes,
                dt_g,
                use_linear,
                &delta_remote,
                plan,
                accum,
                &mut output_integrals,
                &mut input_integrals,
            )?;

            // Advance working accepted states for nodes in this group.
            for &n in &g_nodes {
                core::mem::swap(&mut self.accepted[n], &mut self.candidate[n]);
            }
            t_sub.add(dt_g);

            // Check if a newly reached component event falls before macro end.
            if k + 1 < substeps {
                let elapsed_so_far = ((k + 1) as f64) * dt_g;
                let remaining = macro_dt_s - elapsed_so_far;
                for &n in &g_nodes {
                    if let Some(ev) = self.nodes[n].next_event_s(self.accepted[n].as_ref()) {
                        if !(ev > 0.0) {
                            return Err(SystemError::Negotiation {
                                node: self.names[n].clone(),
                                proposed_s: remaining,
                                returned_s: ev,
                            }
                            .into());
                        }
                        if ev < remaining * (1.0 - 1.0e-12) {
                            return Ok(GroupOutcome::Shorten {
                                macro_dt_s: elapsed_so_far + ev,
                            });
                        }
                    }
                }
            }
        }

        // Compute exact double-double window averages for all outputs/inputs of this group.
        let inv_h = Dd::from_f64(macro_dt_s);
        for (o, integral) in output_integrals {
            accum.window_outputs[o] = (integral / inv_h).to_f64();
        }
        for (i, integral) in input_integrals {
            accum.window_inputs[i] = (integral / inv_h).to_f64();
        }
        Ok(GroupOutcome::Done)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn audit_admit_and_accumulate_substep(
        &mut self,
        g_nodes: &[usize],
        dt_g: f64,
        use_linear: bool,
        delta_remote: &[f64],
        plan: &MultiratePlan,
        accum: &mut MacroAccum,
        output_integrals: &mut BTreeMap<usize, Dd>,
        input_integrals: &mut BTreeMap<usize, Dd>,
    ) -> Result<(), MultirateError> {
        let dt_dd = Dd::from_f64(dt_g);
        let mut delta_u_buf = Vec::new();
        let mut avg_u_buf = Vec::new();
        for &node in g_nodes {
            let stored_old = self.measure_node_energy(node, self.accepted[node].as_ref())?;
            let stored_new = self.measure_node_energy(node, self.candidate[node].as_ref())?;
            let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
            let ports = self.nodes[node].ports();
            let audit = audit_step(
                ports,
                &self.layouts[node],
                dt_g,
                &self.inputs[inputs.clone()],
                &self.outputs[outputs.clone()],
                &self.receipts[node],
                stored_old,
                stored_new,
            )
            .map_err(|violation| SystemError::Contract {
                node: self.names[node].clone(),
                violation,
            })?;

            let mut node_port_e = Dd::ZERO;
            for (port, spec) in ports.iter().enumerate() {
                let in_r = self.layouts[node].input_range(port);
                let out_r = self.layouts[node].output_range(port);
                let input = &self.inputs[inputs.start + in_r.start..inputs.start + in_r.end];
                let output = &self.outputs[outputs.start + out_r.start..outputs.start + out_r.end];
                let global = self.port_start[node] + port;

                let e_j = spec.energy_in_j(dt_g, input, output);
                let s_j_k = spec.entropy_in_j_per_k(dt_g, input, output);
                let m_kg = spec.mass_in_kg(dt_g, input, output);
                node_port_e = node_port_e + Dd::from_f64(e_j);
                accum.port_energy[global] = accum.port_energy[global] + Dd::from_f64(e_j);
                accum.port_entropy[global] = accum.port_entropy[global] + Dd::from_f64(s_j_k);
                accum.port_mass[global] = accum.port_mass[global] + Dd::from_f64(m_kg);
                for (c, slot) in accum.port_carried[global].iter_mut().enumerate() {
                    let amt = spec.carried_quantity_in(dt_g, input, output, c);
                    *slot = *slot + Dd::from_f64(amt);
                }

                let mag = if spec.field_metric().is_some() {
                    dt_g * power_pairing_scale(spec, input, output)
                } else if spec.carrier() == Carrier::Heat && spec.width() > 1 {
                    let loads = if spec.role() == PortRole::EffortIn {
                        output
                    } else {
                        input
                    };
                    dt_g * loads.iter().map(|v| v.abs()).sum::<f64>()
                } else {
                    e_j.abs()
                };
                accum.port_magnitude[global] = accum.port_magnitude[global] + Dd::from_f64(mag);

                // Accumulate solve allowances across substeps.
                for (offset, idx) in
                    (inputs.start + in_r.start..inputs.start + in_r.end).enumerate()
                {
                    let allow = self.input_allowance[idx];
                    if allow > 0.0 {
                        if spec.carrier() == Carrier::Stream
                            && spec.role() == PortRole::FlowIn
                            && offset == STREAM_MASS_RATE
                        {
                            accum.port_mass_slack_kg[global] += dt_g * allow;
                        }
                        if spec.carrier() == Carrier::Stream
                            && spec.role() == PortRole::FlowIn
                            && offset >= STREAM_CARRIED
                        {
                            let c = offset - STREAM_CARRIED;
                            if c < accum.port_carried_slack[global].len() {
                                accum.port_carried_slack[global][c] += dt_g * allow;
                            }
                        }
                        let slack_j = match (spec.carrier(), spec.role()) {
                            (Carrier::Heat, PortRole::FlowIn) => dt_g * allow,
                            (Carrier::Stream, PortRole::FlowIn)
                                if offset == STREAM_ENTHALPY_RATE =>
                            {
                                dt_g * allow
                            }
                            (Carrier::Power(_), _) => {
                                power_input_slack_j(spec, output, offset, dt_g, allow)
                            }
                            _ => 0.0,
                        };
                        accum.port_solve_slack_j[global] += slack_j;
                        if spec.carrier() == Carrier::Heat && spec.role() == PortRole::FlowIn {
                            accum.port_entropy_slack[global] += slack_j / output[offset];
                        }
                    }
                }

                // Exact interpolation quadrature discrepancy for Linear slow-to-fast mode.
                if use_linear {
                    let port_in_range = inputs.start + in_r.start..inputs.start + in_r.end;
                    let has_slower = port_in_range
                        .clone()
                        .any(|idx| plan.input_has_slower_remote[idx]);
                    if has_slower {
                        delta_u_buf.clear();
                        avg_u_buf.clear();
                        for (offset, idx) in port_in_range.enumerate() {
                            let du = driver_remote_delta(
                                &self.drivers[idx],
                                &self.outputs,
                                delta_remote,
                                &accum.window_outputs,
                                &self.boundary_values,
                            );
                            delta_u_buf.push(du);
                            avg_u_buf.push(input[offset] - du);
                        }
                        let interp_e = match spec.carrier() {
                            Carrier::Power(_) => spec.energy_in_j(dt_g, &delta_u_buf, output),
                            _ => 0.0,
                        };
                        if interp_e != 0.0 {
                            let d = Dd::from_f64(interp_e);
                            accum.port_interp_energy[global] = accum.port_interp_energy[global] + d;
                            accum.nodes[node].interp_energy_j =
                                accum.nodes[node].interp_energy_j + d;
                        }
                        if spec.carrier() == Carrier::Heat && spec.role() == PortRole::EffortIn {
                            let s_avg = spec.entropy_in_j_per_k(dt_g, &avg_u_buf, output);
                            let ds = s_j_k - s_avg;
                            if ds != 0.0 {
                                accum.port_interp_entropy[global] =
                                    accum.port_interp_entropy[global] + Dd::from_f64(ds);
                            }
                        }
                    }
                }
            }

            let receipt = &self.receipts[node];
            let acc = &mut accum.nodes[node];
            acc.port_energy_in_j = acc.port_energy_in_j + node_port_e;
            acc.stored_change_j = acc.stored_change_j + Dd::from_f64(receipt.stored_change_j);
            acc.dissipated_j = acc.dissipated_j + Dd::from_f64(receipt.dissipated_j);
            acc.numerical_loss_j = acc.numerical_loss_j + Dd::from_f64(receipt.numerical_loss_j);
            acc.entropy_generated_j_per_k =
                acc.entropy_generated_j_per_k + Dd::from_f64(receipt.entropy_generated_j_per_k);
            acc.balance_defect_j = acc.balance_defect_j + Dd::from_f64(audit.balance_defect_j);
            acc.absolute_defect_j =
                acc.absolute_defect_j + Dd::from_f64(audit.balance_defect_j.abs());

            self.accumulate_node_inventories(node, dt_g, accum)?;

            self.nodes[node]
                .admit(
                    self.accepted[node].as_ref(),
                    self.candidate[node].as_ref(),
                    &self.receipts[node],
                )
                .map_err(|error| SystemError::Component {
                    node: self.names[node].clone(),
                    phase: Phase::Admit,
                    error,
                })?;

            for o in self.node_outputs(node) {
                if let Some(slot) = output_integrals.get_mut(&o) {
                    *slot = *slot + Dd::from_f64(self.outputs[o]) * dt_dd;
                }
            }
            for i in self.node_inputs(node) {
                if let Some(slot) = input_integrals.get_mut(&i) {
                    *slot = *slot + Dd::from_f64(self.inputs[i]) * dt_dd;
                }
                accum.window_allowance[i] = accum.window_allowance[i].max(self.input_allowance[i]);
            }
        }
        Ok(())
    }

    fn accumulate_node_inventories(
        &self,
        node: usize,
        dt_g: f64,
        accum: &mut MacroAccum,
    ) -> Result<(), MultirateError> {
        let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
        let ports = self.nodes[node].ports();

        if self.ledger.nodes[node].has_streams {
            let old_m = self.measure_node_mass(node, self.accepted[node].as_ref())?;
            let new_m = self.measure_node_mass(node, self.candidate[node].as_ref())?;
            let mass_audit: Option<MassStepAudit> = audit_mass_step(
                ports,
                &self.layouts[node],
                dt_g,
                &self.inputs[inputs.clone()],
                &self.outputs[outputs.clone()],
                old_m,
                new_m,
            )
            .map_err(|violation| SystemError::Contract {
                node: self.names[node].clone(),
                violation,
            })?;
            if let Some(m) = mass_audit {
                let acc = &mut accum.nodes[node];
                acc.has_mass = true;
                acc.mass_port_in_kg = acc.mass_port_in_kg + Dd::from_f64(m.port_mass_in_kg);
                acc.mass_stored_change_kg =
                    acc.mass_stored_change_kg + Dd::from_f64(m.stored_change_kg);
                acc.mass_balance_defect_kg =
                    acc.mass_balance_defect_kg + Dd::from_f64(m.balance_defect_kg);
                acc.mass_abs_defect_kg =
                    acc.mass_abs_defect_kg + Dd::from_f64(m.balance_defect_kg.abs());
                acc.mass_audited_steps += 1;
            }
        }

        let entropy_audit: Option<EntropyStepAudit> = self.audit_entropy(node, dt_g)?;
        if let Some(s) = entropy_audit {
            let acc = &mut accum.nodes[node];
            acc.has_entropy = true;
            acc.entropy_port_in = acc.entropy_port_in + Dd::from_f64(s.port_entropy_in_j_per_k);
            acc.entropy_gen_audit =
                acc.entropy_gen_audit + Dd::from_f64(s.entropy_generated_j_per_k);
            acc.entropy_stored_change =
                acc.entropy_stored_change + Dd::from_f64(s.stored_change_j_per_k);
            acc.entropy_balance_defect =
                acc.entropy_balance_defect + Dd::from_f64(s.balance_defect_j_per_k);
            acc.entropy_abs_defect =
                acc.entropy_abs_defect + Dd::from_f64(s.balance_defect_j_per_k.abs());
            acc.entropy_audited_steps += 1;
        }

        for (c, slot_opt) in accum.nodes[node].carried.iter_mut().enumerate() {
            let old_c = self.nodes[node]
                .stored_carried_quantity(self.accepted[node].as_ref(), c)
                .map_err(|error| SystemError::Component {
                    node: self.names[node].clone(),
                    phase: Phase::Measure,
                    error,
                })?;
            let new_c = self.nodes[node]
                .stored_carried_quantity(self.candidate[node].as_ref(), c)
                .map_err(|error| SystemError::Component {
                    node: self.names[node].clone(),
                    phase: Phase::Measure,
                    error,
                })?;
            let ca: Option<CarriedStepAudit> = audit_carried_step(
                ports,
                &self.layouts[node],
                dt_g,
                &self.inputs[inputs.clone()],
                &self.outputs[outputs.clone()],
                c,
                old_c,
                new_c,
            )
            .map_err(|violation| SystemError::Contract {
                node: self.names[node].clone(),
                violation,
            })?;
            if let Some(a) = ca {
                let slot = slot_opt.get_or_insert_with(Default::default);
                slot.port_in = slot.port_in + Dd::from_f64(a.port_in);
                slot.stored_change = slot.stored_change + Dd::from_f64(a.stored_change);
                slot.balance_defect = slot.balance_defect + Dd::from_f64(a.balance_defect);
                slot.abs_defect = slot.abs_defect + Dd::from_f64(a.balance_defect.abs());
                slot.audited_steps += 1;
            }
        }
        Ok(())
    }

    fn measure_node_energy(&self, node: usize, state: &dyn Any) -> Result<f64, SystemError> {
        self.nodes[node]
            .stored_energy_j(state)
            .map_err(|error| SystemError::Component {
                node: self.names[node].clone(),
                phase: Phase::Measure,
                error,
            })
    }

    fn measure_node_mass(&self, node: usize, state: &dyn Any) -> Result<Option<f64>, SystemError> {
        self.nodes[node]
            .stored_mass_kg(state)
            .map_err(|error| SystemError::Component {
                node: self.names[node].clone(),
                phase: Phase::Measure,
                error,
            })
    }

    #[allow(clippy::too_many_lines)]
    fn solve_local_sub_block(
        &mut self,
        sub: &LocalSubBlock,
        dt_g: f64,
        cx: &Cx<'_>,
    ) -> Result<SubBlockOutcome, MultirateError> {
        let m = sub.unknowns.len();
        for &idx in &sub.unknowns {
            if self.input_is_temperature[idx] && !(self.inputs[idx] > 0.0) {
                return Err(SystemError::GuessRequired {
                    port: self.input_label[idx].clone(),
                }
                .into());
            }
        }

        let mut condensed = Condensed::plan(&sub.condensed_pattern);
        let mut residual = vec![0.0_f64; m];
        let mut group_scale = vec![0.0_f64; self.group_count];
        let mut rounding = vec![0.0_f64; m];
        let mut allowance = vec![0.0_f64; m];
        let mut grad = Vec::new();
        let mut iterations = 0_u32;

        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            for &node in &sub.nodes {
                if let Some(shorter) = self.propose(node, dt_g, cx)? {
                    return Ok(SubBlockOutcome::Shorten { dt_s: shorter });
                }
            }
            group_scale.fill(0.0);
            for (row, &idx) in sub.unknowns.iter().enumerate() {
                let drv = &self.drivers[idx];
                let target = drv.value(&self.outputs, &self.boundary_values);
                residual[row] = self.inputs[idx] - target;
                let g = self.input_group[idx];
                let own = self.inputs[idx].abs();
                let terms = drv.magnitude(&self.outputs, &self.boundary_values);
                group_scale[g] = group_scale[g].max(own).max(terms);
            }

            let power_scale = self.block_power_scale(&sub.nodes);
            let mut converged = true;
            let mut worst_row = 0_usize;
            let mut worst_ratio = -1.0_f64;
            for (row, &idx) in sub.unknowns.iter().enumerate() {
                let mut tol_allow =
                    self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * rounding[row];
                if self.input_is_energy_rate[idx] {
                    tol_allow += ROUNDING_UNITS * f64::EPSILON * power_scale;
                }
                allowance[row] = tol_allow;
                let g = self.input_group[idx];
                let tol = self.config.relative_tolerance * group_scale[g] + tol_allow;
                if !(residual[row].is_finite() && residual[row].abs() <= tol) {
                    converged = false;
                    let ratio = if tol > 0.0 {
                        residual[row].abs() / tol
                    } else {
                        f64::INFINITY
                    };
                    if !ratio.is_finite() || ratio > worst_ratio {
                        worst_ratio = ratio;
                        worst_row = row;
                    }
                }
            }
            if converged {
                for (row, &idx) in sub.unknowns.iter().enumerate() {
                    self.input_allowance[idx] = allowance[row];
                }
                return Ok(SubBlockOutcome::Done { iterations });
            }
            if iterations >= self.config.max_iterations {
                return Err(self
                    .not_converged(&sub.nodes, &sub.unknowns, &residual, iterations, worst_row)
                    .into());
            }

            // Assemble condensed Jacobian I - d(driver)/d(unknown).
            condensed.clear();
            rounding.fill(0.0);
            let mut row_grads: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
            for &idx in &sub.unknowns {
                self.drivers[idx].gradient(&self.outputs, &mut grad);
                row_grads.push(grad.clone());
            }
            for (col, &idx) in sub.unknowns.iter().enumerate() {
                let owner = self.input_owner[idx];
                let local = idx - self.input_start[owner];
                if let Some(shorter) =
                    self.tangent_column(owner, local, dt_g, &group_scale, false, cx)?
                {
                    return Ok(SubBlockOutcome::Shorten { dt_s: shorter });
                }
                let input_mag = self.inputs[idx].abs();
                for (local_out, out_idx) in self.node_outputs(owner).enumerate() {
                    let dy = self.scratch_outputs[local_out];
                    if dy == 0.0 {
                        continue;
                    }
                    if let Some(rows) = sub.readers.get(&out_idx) {
                        for &row in rows {
                            let coeff = row_grads[row]
                                .iter()
                                .find(|&&(o, _)| o == out_idx)
                                .map_or(0.0, |&(_, c)| c);
                            let entry = coeff * dy;
                            if entry != 0.0 {
                                let slot = condensed.slot(row, col);
                                condensed.subtract(slot, entry);
                                rounding[row] += entry.abs() * input_mag;
                            }
                        }
                    }
                }
            }
            if !condensed.factor() {
                return Err(SystemError::Singular {
                    nodes: self.names_of(&sub.nodes),
                }
                .into());
            }
            for v in &mut residual {
                *v = -*v;
            }
            if !condensed.solve(&mut residual) {
                return Err(SystemError::Singular {
                    nodes: self.names_of(&sub.nodes),
                }
                .into());
            }

            // Apply damped update preserving strict positivity on temperature inputs.
            let mut alpha = 1.0_f64;
            for _ in 0..16 {
                let temp_ok = sub.unknowns.iter().enumerate().all(|(r, &idx)| {
                    !self.input_is_temperature[idx]
                        || (self.inputs[idx] + alpha * residual[r]) > 0.0
                });
                if temp_ok {
                    break;
                }
                alpha *= 0.5;
            }
            for (row, &idx) in sub.unknowns.iter().enumerate() {
                self.inputs[idx] += alpha * residual[row];
            }
            iterations += 1;
        }
    }

    fn solve_cyclic_stage(
        &mut self,
        stage: &GroupStage,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        entry: &Checkpoint,
        accum: &mut MacroAccum,
        cx: &Cx<'_>,
    ) -> Result<CyclicStageOutcome, MultirateError> {
        let k = stage.back_outputs.len();
        let mut z = self.seed_cyclic_back_outputs(stage, macro_dt_s, plan, entry, accum, cx)?;
        let mut r = vec![0.0_f64; k];
        let mut y_new = vec![0.0_f64; k];
        let mut jac = vec![0.0_f64; k * k];

        for iter in 0..=plan.config.cross_rate_max_iterations {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            match self.eval_cyclic_stage_once(stage, &z, macro_dt_s, plan, entry, accum, cx)? {
                GroupOutcome::Shorten { macro_dt_s: s } => {
                    return Ok(CyclicStageOutcome::Shorten { macro_dt_s: s });
                }
                GroupOutcome::Done => {}
            }
            let mut converged = true;
            let mut worst_rel = 0.0_f64;
            for (i, &out_idx) in stage.back_outputs.iter().enumerate() {
                y_new[i] = accum.window_outputs[out_idx];
                r[i] = z[i] - y_new[i];
                let scale = z[i].abs().max(y_new[i].abs()).max(1.0e-30);
                let tol = plan.config.cross_rate_relative_tolerance * scale
                    + plan.config.cross_rate_absolute_tolerance
                    + ROUNDING_UNITS * f64::EPSILON * scale;
                let rel = r[i].abs() / scale;
                if rel > worst_rel {
                    worst_rel = rel;
                }
                if !(r[i].is_finite() && r[i].abs() <= tol) {
                    converged = false;
                }
            }
            if converged {
                return Ok(CyclicStageOutcome::Done {
                    iterations: iter,
                    max_residual: worst_rel,
                });
            }
            if iter == plan.config.cross_rate_max_iterations {
                return Err(MultirateError::CrossRateNotConverged {
                    iterations: iter,
                    residual: worst_rel,
                });
            }

            // Finite-difference Jacobian of R(z) = z - y_new(z).
            for col in 0..k {
                cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
                let mut z_p = z.clone();
                let h = 1.0e-6 * (1.0 + z[col].abs());
                z_p[col] += h;
                match self
                    .eval_cyclic_stage_once(stage, &z_p, macro_dt_s, plan, entry, accum, cx)?
                {
                    GroupOutcome::Shorten { macro_dt_s: s } => {
                        return Ok(CyclicStageOutcome::Shorten { macro_dt_s: s });
                    }
                    GroupOutcome::Done => {}
                }
                for (row, &out_idx) in stage.back_outputs.iter().enumerate() {
                    let r_p = z_p[row] - accum.window_outputs[out_idx];
                    let entry_val = (r_p - r[row]) / h;
                    if !entry_val.is_finite() {
                        return Err(MultirateError::CrossRateSingular);
                    }
                    jac[row * k + col] = entry_val;
                }
            }
            let factor = lu(&jac, k).map_err(|_| MultirateError::CrossRateSingular)?;
            let mut dz: Vec<f64> = r.iter().map(|&v| -v).collect();
            factor.solve(&mut dz);
            if dz.iter().any(|v| !v.is_finite()) {
                return Err(MultirateError::CrossRateSingular);
            }
            let mut alpha = 1.0_f64;
            for _ in 0..16 {
                let ok = stage.back_outputs.iter().enumerate().all(|(i, &o)| {
                    !plan.output_feeds_temperature[o] || (z[i] + alpha * dz[i]) > 0.0
                });
                if ok {
                    break;
                }
                alpha *= 0.5;
            }
            for i in 0..k {
                z[i] += alpha * dz[i];
            }
        }
        Err(MultirateError::CrossRateSingular)
    }

    fn seed_cyclic_back_outputs(
        &mut self,
        stage: &GroupStage,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        entry: &Checkpoint,
        accum: &mut MacroAccum,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, MultirateError> {
        let mut z = Vec::with_capacity(stage.back_outputs.len());
        let mut need_slow_probe = false;
        for &o in &stage.back_outputs {
            let mut init = entry.outputs[o];
            if plan.output_feeds_temperature[o] && !(init > 0.0) {
                // Check whether any reader input received an explicit guess.
                let mut found = None;
                let mut src_buf = Vec::new();
                for (in_idx, drv) in self.drivers.iter().enumerate() {
                    if self.input_is_temperature[in_idx] && entry.inputs[in_idx] > 0.0 {
                        src_buf.clear();
                        drv.sources(&mut src_buf);
                        if src_buf.contains(&o) {
                            found = Some(entry.inputs[in_idx]);
                            break;
                        }
                    }
                }
                if let Some(val) = found {
                    init = val;
                } else {
                    need_slow_probe = true;
                }
            }
            z.push(init);
        }

        if need_slow_probe {
            // Probe slower groups once at zero cross-rate input to read their
            // resting outputs (e.g., initial temperature of a HeatCapacity).
            for &g in stage.groups.iter().skip(1) {
                let _ = self.eval_rate_group(g, macro_dt_s, plan, entry, accum, cx);
                let g_nodes: Vec<usize> = plan.groups[g].nodes.iter().map(|id| id.0).collect();
                self.restore_group_nodes(&g_nodes, entry)?;
            }
            for (i, &o) in stage.back_outputs.iter().enumerate() {
                if plan.output_feeds_temperature[o]
                    && !(z[i] > 0.0)
                    && accum.window_outputs[o] > 0.0
                {
                    z[i] = accum.window_outputs[o];
                }
            }
        }
        Ok(z)
    }

    #[allow(clippy::too_many_arguments)]
    fn eval_cyclic_stage_once(
        &mut self,
        stage: &GroupStage,
        z: &[f64],
        macro_dt_s: f64,
        plan: &MultiratePlan,
        entry: &Checkpoint,
        accum: &mut MacroAccum,
        cx: &Cx<'_>,
    ) -> Result<GroupOutcome, MultirateError> {
        for (&o, &val) in stage.back_outputs.iter().zip(z) {
            accum.window_outputs[o] = val;
        }
        for &g in &stage.groups {
            match self.eval_rate_group(g, macro_dt_s, plan, entry, accum, cx)? {
                GroupOutcome::Done => {}
                shorten @ GroupOutcome::Shorten { .. } => return Ok(shorten),
            }
        }
        Ok(GroupOutcome::Done)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn audit_and_publish_macro(
        &mut self,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        accum: &MacroAccum,
        cross_rate_iterations: u32,
        shortenings: u32,
        max_cross_rate_residual: f64,
        cx: &Cx<'_>,
    ) -> Result<MultirateStepReport, MultirateError> {
        let mut members = Vec::new();
        let mut group_members = Vec::new();
        let mut connection_rows = Vec::with_capacity(self.audit_groups.len());
        let mut sink_heat = vec![0.0_f64; self.ledger.sinks.len()];
        let mut generated = Dd::ZERO;
        let mut entropy_scale = 0.0_f64;
        let mut entropy_slack = 0.0_f64;
        let mut total_interp_energy = Dd::ZERO;
        let mut total_interp_entropy = Dd::ZERO;

        for node_acc in &accum.nodes {
            let eg = node_acc.entropy_generated_j_per_k.to_f64();
            generated = generated + node_acc.entropy_generated_j_per_k;
            entropy_scale += eg.abs();
        }

        for (index, group) in self.audit_groups.iter().enumerate() {
            if group.is_empty() {
                connection_rows.push((0.0, 0.0, 0.0));
                self.connection_carried[index].fill(0.0);
                continue;
            }
            members.clear();
            for &connection in group {
                self.connections[connection].energy_members(&mut group_members);
                for port in &group_members {
                    if !members.contains(port) {
                        members.push(*port);
                    }
                }
            }
            let mut energy_dd = Dd::ZERO;
            let mut interp_e_dd = Dd::ZERO;
            let mut interp_s_dd = Dd::ZERO;
            let mut magnitude = 0.0_f64;
            let mut entropy_dd = Dd::ZERO;
            let mut mass_dd = Dd::ZERO;
            let mut mass_magnitude = 0.0_f64;
            let mut solve_slack_j = 0.0_f64;
            let mut mass_slack_kg = 0.0_f64;

            for port in &members {
                let global = self.port_start[port.node] + port.port;
                energy_dd = energy_dd + accum.port_energy[global];
                interp_e_dd = interp_e_dd + accum.port_interp_energy[global];
                interp_s_dd = interp_s_dd + accum.port_interp_entropy[global];
                magnitude += accum.port_magnitude[global].to_f64();
                entropy_dd = entropy_dd + accum.port_entropy[global];
                entropy_scale += accum.port_entropy[global].to_f64().abs();
                mass_dd = mass_dd + accum.port_mass[global];
                mass_magnitude += accum.port_mass[global].to_f64().abs();
                solve_slack_j += accum.port_solve_slack_j[global];
                mass_slack_kg += accum.port_mass_slack_kg[global];
                entropy_slack += accum.port_entropy_slack[global];
            }

            for &connection in group {
                if let Connection::HeatRoute {
                    source,
                    remainder: Some((sink, fraction, kelvin)),
                    ..
                } = &self.connections[connection]
                {
                    let emitted = -accum.port_energy[self.port_start[source.node] + source.port];
                    let heat_dd = Dd::from_f64(*fraction) * emitted;
                    let heat = heat_dd.to_f64();
                    sink_heat[*sink] += heat;
                    energy_dd = energy_dd + heat_dd;
                    magnitude += heat.abs();
                    let s_dd = heat_dd / Dd::from_f64(*kelvin);
                    entropy_dd = entropy_dd + s_dd;
                    entropy_scale += s_dd.to_f64().abs();
                }
            }

            // Separate true conservation defect from explicitly booked
            // slow-to-fast linear interpolation quadrature discrepancy.
            let true_energy = (energy_dd - interp_e_dd).to_f64();
            total_interp_energy = total_interp_energy + interp_e_dd;
            total_interp_entropy = total_interp_entropy + interp_s_dd;

            let tolerance_j = self.config.interface_energy_tolerance * magnitude
                + solve_slack_j
                + f64::MIN_POSITIVE;
            if !(true_energy.is_finite()
                && tolerance_j.is_finite()
                && true_energy.abs() <= tolerance_j)
            {
                return Err(SystemError::Interface {
                    connection: index,
                    defect_j: true_energy,
                    tolerance_j,
                }
                .into());
            }

            let true_entropy = (entropy_dd - interp_s_dd).to_f64();
            generated = generated + (entropy_dd - interp_s_dd);

            let mass_defect_kg = mass_dd.to_f64();
            let mass_tolerance_kg = (128.0 * f64::EPSILON + self.config.relative_tolerance)
                * mass_magnitude
                + mass_slack_kg;
            if !(mass_defect_kg.is_finite()
                && mass_tolerance_kg.is_finite()
                && mass_defect_kg.abs() <= mass_tolerance_kg)
            {
                return Err(SystemError::InterfaceMass {
                    connection: index,
                    defect_kg: mass_defect_kg,
                    tolerance_kg: mass_tolerance_kg,
                }
                .into());
            }
            connection_rows.push((true_energy, true_entropy, mass_defect_kg));

            for carried in 0..self.connection_carried[index].len() {
                let mut c_sum = Dd::ZERO;
                let mut c_mag = 0.0_f64;
                let mut c_slack = 0.0_f64;
                for port in &members {
                    let global = self.port_start[port.node] + port.port;
                    if let Some(&amt) = accum.port_carried[global].get(carried) {
                        c_sum = c_sum + amt;
                        c_mag += amt.to_f64().abs();
                    }
                    if let Some(&sl) = accum.port_carried_slack[global].get(carried) {
                        c_slack += sl;
                    }
                }
                let defect = c_sum.to_f64();
                let tolerance =
                    (128.0 * f64::EPSILON + self.config.relative_tolerance) * c_mag + c_slack;
                if !(defect.is_finite() && tolerance.is_finite() && defect.abs() <= tolerance) {
                    return Err(SystemError::InterfaceCarried {
                        connection: index,
                        carried,
                        defect,
                        tolerance,
                    }
                    .into());
                }
                self.connection_carried[index][carried] = defect;
            }
        }

        let generated_f64 = generated.to_f64();
        if !(generated_f64 >= -1.0e-9 * entropy_scale - entropy_slack) {
            return Err(SystemError::SecondLaw {
                generated_j_per_k: generated_f64,
            }
            .into());
        }

        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;

        // Atomic publication: all states are already in `self.accepted`.
        self.accepted_inputs.copy_from_slice(&accum.window_inputs);
        self.accepted_outputs.copy_from_slice(&accum.window_outputs);
        self.accepted_allowance
            .copy_from_slice(&accum.window_allowance);
        self.accepted_receipts.copy_from_slice(&self.receipts);
        self.inputs.copy_from_slice(&self.accepted_inputs);
        self.outputs.copy_from_slice(&self.accepted_outputs);
        self.accepted_generation = self.accepted_generation.wrapping_add(1);
        self.invalidate_caches();

        self.time.add(macro_dt_s);
        self.steps += 1;
        self.ledger.steps += 1;
        self.ledger.duration_s.add(macro_dt_s);

        for (node, account) in self.ledger.nodes.iter_mut().enumerate() {
            let acc = &accum.nodes[node];
            let g_substeps = plan.groups[plan.node_group[node]].substeps as u64;
            // Book exact slow-to-fast interpolation discrepancy as numerical
            // loss while keeping both per-node and global energy identities exact.
            let net_port_e = (acc.port_energy_in_j - acc.interp_energy_j).to_f64();
            let net_num_loss = (acc.numerical_loss_j - acc.interp_energy_j).to_f64();
            account.port_energy_in_j.add(net_port_e);
            account.stored_change_j.add(acc.stored_change_j.to_f64());
            account.dissipated_j.add(acc.dissipated_j.to_f64());
            account.numerical_loss_j.add(net_num_loss);
            account
                .entropy_generated_j_per_k
                .add(acc.entropy_generated_j_per_k.to_f64());
            account.balance_defect_j.add(acc.balance_defect_j.to_f64());
            account
                .absolute_defect_j
                .add(acc.absolute_defect_j.to_f64());
            if acc.has_mass {
                let mass = account.mass.get_or_insert_with(Default::default);
                mass.port_mass_in_kg.add(acc.mass_port_in_kg.to_f64());
                mass.stored_change_kg
                    .add(acc.mass_stored_change_kg.to_f64());
                mass.balance_defect_kg
                    .add(acc.mass_balance_defect_kg.to_f64());
                mass.absolute_defect_kg.add(acc.mass_abs_defect_kg.to_f64());
                if acc.mass_audited_steps == g_substeps {
                    mass.audited_steps += 1;
                }
            }
            if acc.has_entropy {
                let ent = account.entropy.get_or_insert_with(Default::default);
                ent.port_entropy_in_j_per_k
                    .add(acc.entropy_port_in.to_f64());
                ent.entropy_generated_j_per_k
                    .add(acc.entropy_gen_audit.to_f64());
                ent.stored_change_j_per_k
                    .add(acc.entropy_stored_change.to_f64());
                ent.balance_defect_j_per_k
                    .add(acc.entropy_balance_defect.to_f64());
                ent.absolute_defect_j_per_k
                    .add(acc.entropy_abs_defect.to_f64());
                if acc.entropy_audited_steps == g_substeps {
                    ent.audited_steps += 1;
                }
            }
            for (slot, maybe_ca) in account.carried.iter_mut().zip(&acc.carried) {
                if let Some(ca) = maybe_ca {
                    let target = slot.get_or_insert_with(Default::default);
                    target.port_in.add(ca.port_in.to_f64());
                    target.stored_change.add(ca.stored_change.to_f64());
                    target.balance_defect.add(ca.balance_defect.to_f64());
                    target.absolute_defect.add(ca.abs_defect.to_f64());
                    if ca.audited_steps == g_substeps {
                        target.audited_steps += 1;
                    }
                }
            }
        }
        for (account, &(energy, entropy, mass)) in
            self.ledger.connections.iter_mut().zip(&connection_rows)
        {
            account.energy_defect_j.add(energy);
            account.entropy_generated_j_per_k.add(entropy);
            account.mass_defect_kg.add(mass);
            account.absolute_energy_defect_j.add(energy.abs());
            account.absolute_mass_defect_kg.add(mass.abs());
        }
        for (account, row) in self
            .ledger
            .connections
            .iter_mut()
            .zip(&self.connection_carried)
        {
            for (total, &value) in account.carried_defect.iter_mut().zip(row) {
                total.add(value);
            }
        }
        for (account, &heat) in self.ledger.sinks.iter_mut().zip(&sink_heat) {
            account.heat_j.add(heat);
        }
        for (account, (port, _)) in self.ledger.boundaries.iter_mut().zip(&self.boundaries) {
            let global = self.port_start[port.node] + port.port;
            account.energy_in_j.add(accum.port_energy[global].to_f64());
            account
                .entropy_in_j_per_k
                .add(accum.port_entropy[global].to_f64());
            account.mass_in_kg.add(accum.port_mass[global].to_f64());
            for (total, value) in account
                .carried_in
                .iter_mut()
                .zip(&accum.port_carried[global])
            {
                total.add(value.to_f64());
            }
        }

        let group_substeps = plan.groups.iter().map(|g| g.substeps).collect();
        Ok(MultirateStepReport {
            dt_s: macro_dt_s,
            time_s: self.time.value(),
            group_substeps,
            cross_rate_iterations,
            inner_iterations: accum.inner_iterations,
            shortenings,
            interpolation_numerical_loss_j: (-total_interp_energy).to_f64(),
            interpolation_entropy_defect_j_per_k: (-total_interp_entropy).to_f64(),
            max_cross_rate_residual,
        })
    }

    /// Advance one multirate macro step to the earliest requested zero crossing
    /// among `guards`, synchronizing every rate group to the exact same clock.
    ///
    /// When `plan.is_single_rate()` is `true`, delegates directly to
    /// [`Self::step_to_events`] for bitwise identity.
    ///
    /// # Errors
    /// Foreign plan, invalid event config or guard output, failed/cancelled
    /// multirate step, or exhausted/unresolvable bisection bracket. Every error
    /// restores the full entry state and accounts.
    #[allow(clippy::too_many_lines)]
    pub fn step_multirate_to_events(
        &mut self,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        config: EventSetConfig,
        directions: &[EventDirection],
        mut guards: impl FnMut(&System, &mut [f64]),
        cx: &Cx<'_>,
    ) -> Result<MultirateEventReport, MultirateError> {
        if !Rc::ptr_eq(&self.identity, &plan.identity) {
            return Err(MultirateError::ForeignPlan);
        }
        if plan.is_single_rate() {
            let report = self.step_to_events(macro_dt_s, config, directions, guards, cx)?;
            let step = report.step.map(|s: StepReport| MultirateStepReport {
                dt_s: s.dt_s,
                time_s: s.time_s,
                group_substeps: vec![1; plan.groups.len()],
                cross_rate_iterations: 0,
                inner_iterations: s.iterations,
                shortenings: s.shortenings,
                interpolation_numerical_loss_j: 0.0,
                interpolation_entropy_defect_j_per_k: 0.0,
                max_cross_rate_residual: 0.0,
            });
            return Ok(MultirateEventReport {
                status: report.status,
                step,
                time_s: report.time_s,
                values: report.values,
                trial_steps: report.trial_steps,
            });
        }
        if directions.is_empty() {
            return Err(EventSetError::EmptySet.into());
        }
        if !(config.time_tolerance_s.is_finite() && config.time_tolerance_s > 0.0) {
            return Err(EventSetError::InvalidConfig {
                field: "time_tolerance_s",
            }
            .into());
        }
        if config.max_bisections == 0 {
            return Err(EventSetError::InvalidConfig {
                field: "max_bisections",
            }
            .into());
        }
        if !(macro_dt_s.is_finite() && macro_dt_s > 0.0) {
            return Err(SystemError::Step { dt_s: macro_dt_s }.into());
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;

        let entry = self.checkpoint()?;
        let entry_allowance = self.accepted_allowance.clone();
        let entry_receipts = self.accepted_receipts.clone();
        let result = self.multirate_event_trials(
            macro_dt_s,
            plan,
            config,
            directions,
            &entry,
            &entry_allowance,
            &entry_receipts,
            &mut guards,
            cx,
        );
        if result.is_err() {
            self.restore(&entry)?;
            self.accepted_allowance.copy_from_slice(&entry_allowance);
            self.accepted_receipts.copy_from_slice(&entry_receipts);
        }
        result
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn multirate_event_trials(
        &mut self,
        macro_dt_s: f64,
        plan: &MultiratePlan,
        config: EventSetConfig,
        directions: &[EventDirection],
        entry: &Checkpoint,
        entry_allowance: &[f64],
        entry_receipts: &[fs_ports::component::Receipt],
        guards: &mut impl FnMut(&System, &mut [f64]),
        cx: &Cx<'_>,
    ) -> Result<MultirateEventReport, MultirateError> {
        let count = directions.len();
        let mut before = vec![0.0_f64; count];
        eval_guards(self, guards, &mut before, cx)?;
        let zeros: Vec<usize> = before
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| (v == 0.0).then_some(i))
            .collect();
        if !zeros.is_empty() {
            return Ok(MultirateEventReport {
                status: EventSetStatus::AtStart(zeros),
                step: None,
                time_s: self.time_s(),
                values: before,
                trial_steps: 0,
            });
        }

        let mut upper_step = self.step_multirate(macro_dt_s, plan, cx)?;
        let mut upper_dt = upper_step.dt_s;
        let mut upper_values = vec![0.0_f64; count];
        eval_guards(self, guards, &mut upper_values, cx)?;
        let mut trial_steps = 1_u64;

        let crossed_any = directions
            .iter()
            .zip(&before)
            .zip(&upper_values)
            .any(|((dir, &b), &a)| dir_crosses(*dir, b, a));
        if !crossed_any {
            return Ok(MultirateEventReport {
                status: EventSetStatus::Advanced,
                step: Some(upper_step),
                time_s: self.time_s(),
                values: upper_values,
                trial_steps,
            });
        }

        let mut lower_dt = 0.0_f64;
        let mut lower_values = before.clone();
        let mut upper_ckpt = self.checkpoint()?;
        let mut upper_allow = self.accepted_allowance.clone();
        let mut upper_receipts = self.accepted_receipts.clone();
        let mut trial_values = vec![0.0_f64; count];
        let mut bisections = 0_u32;

        while upper_dt - lower_dt > config.time_tolerance_s {
            if bisections >= config.max_bisections {
                let events = collect_located_events(
                    directions,
                    &before,
                    &lower_values,
                    &upper_values,
                    lower_dt,
                    upper_dt,
                );
                return Err(EventSetError::Bisections {
                    count: bisections,
                    events,
                }
                .into());
            }
            let mid_dt = lower_dt + 0.5 * (upper_dt - lower_dt);
            if !(mid_dt > lower_dt && mid_dt < upper_dt) {
                return Err(EventSetError::Resolution {
                    lower_dt_s: lower_dt,
                    upper_dt_s: upper_dt,
                }
                .into());
            }

            self.restore(entry)?;
            self.accepted_allowance.copy_from_slice(entry_allowance);
            self.accepted_receipts.copy_from_slice(entry_receipts);
            let step = self.step_multirate(mid_dt, plan, cx)?;
            trial_steps += 1;
            bisections += 1;
            let actual_dt = step.dt_s;
            if !(actual_dt > lower_dt && actual_dt < upper_dt) {
                return Err(EventSetError::Resolution {
                    lower_dt_s: lower_dt,
                    upper_dt_s: upper_dt,
                }
                .into());
            }
            eval_guards(self, guards, &mut trial_values, cx)?;
            let mid_crossed: Vec<usize> = directions
                .iter()
                .zip(&before)
                .zip(&trial_values)
                .enumerate()
                .filter_map(|(i, ((dir, &b), &a))| dir_crosses(*dir, b, a).then_some(i))
                .collect();
            if mid_crossed.is_empty() {
                lower_dt = actual_dt;
                lower_values.copy_from_slice(&trial_values);
            } else {
                upper_dt = actual_dt;
                upper_values.copy_from_slice(&trial_values);
                upper_step = step;
                upper_ckpt = self.checkpoint()?;
                upper_allow.copy_from_slice(&self.accepted_allowance);
                upper_receipts.copy_from_slice(&self.accepted_receipts);
                if mid_crossed.len() == 1 && upper_values[mid_crossed[0]] == 0.0 {
                    lower_dt = upper_dt;
                    lower_values.copy_from_slice(&upper_values);
                    break;
                }
            }
        }

        self.restore(&upper_ckpt)?;
        self.accepted_allowance.copy_from_slice(&upper_allow);
        self.accepted_receipts.copy_from_slice(&upper_receipts);
        let events = collect_located_events(
            directions,
            &before,
            &lower_values,
            &upper_values,
            lower_dt,
            upper_dt,
        );
        Ok(MultirateEventReport {
            status: EventSetStatus::Located(events),
            step: Some(upper_step),
            time_s: self.time_s(),
            values: upper_values,
            trial_steps,
        })
    }
}

enum CyclicStageOutcome {
    Done { iterations: u32, max_residual: f64 },
    Shorten { macro_dt_s: f64 },
}

fn dir_crosses(dir: EventDirection, before: f64, after: f64) -> bool {
    match dir {
        EventDirection::Either => (before < 0.0 && after >= 0.0) || (before > 0.0 && after <= 0.0),
        EventDirection::Rising => before < 0.0 && after >= 0.0,
        EventDirection::Falling => before > 0.0 && after <= 0.0,
    }
}

fn eval_guards(
    system: &System,
    guards: &mut impl FnMut(&System, &mut [f64]),
    values: &mut [f64],
    cx: &Cx<'_>,
) -> Result<(), MultirateError> {
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    values.fill(f64::NAN);
    guards(system, values);
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    for (index, &value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(EventSetError::InvalidGuard {
                index,
                time_s: system.time_s(),
                value,
            }
            .into());
        }
    }
    Ok(())
}

fn collect_located_events(
    directions: &[EventDirection],
    before: &[f64],
    lower_values: &[f64],
    upper_values: &[f64],
    lower_dt_s: f64,
    upper_dt_s: f64,
) -> Vec<LocatedEvent> {
    directions
        .iter()
        .zip(before)
        .zip(upper_values)
        .enumerate()
        .filter_map(|(index, ((dir, &b), &u))| {
            dir_crosses(*dir, b, u).then_some(LocatedEvent {
                index,
                bracket: EventBracket {
                    lower_dt_s,
                    upper_dt_s,
                    lower_value: lower_values[index],
                    upper_value: u,
                },
            })
        })
        .collect()
}

/// Exact difference `driver(outputs) - driver(outputs_with_remote_at_window_avg)`
/// for one scalar input driven by linearly interpolated slower outputs.
fn driver_remote_delta(
    driver: &Driver,
    outputs: &[f64],
    delta_remote: &[f64],
    window_outputs: &[f64],
    boundaries: &[Vec<f64>],
) -> f64 {
    match driver {
        Driver::Linear(terms) | Driver::OnLine { terms, .. } => {
            let mut sum = 0.0;
            for &(out_idx, coeff) in terms {
                let d = delta_remote[out_idx];
                if d != 0.0 {
                    sum += coeff * d;
                }
            }
            sum
        }
        Driver::EquivalentTemperature {
            terms,
            fixed,
            total,
        } => {
            let cur = driver.value(outputs, boundaries);
            let mut inv_avg = *fixed;
            for &(out_idx, frac) in terms {
                let t = if delta_remote[out_idx] == 0.0 {
                    outputs[out_idx]
                } else {
                    window_outputs[out_idx]
                };
                inv_avg += frac / t;
            }
            cur - (*total / inv_avg)
        }
        Driver::Unset | Driver::Boundary { .. } | Driver::Fixed(_) => 0.0,
    }
}
