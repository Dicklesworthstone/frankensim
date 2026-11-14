//! The system: components, connections, the state store and the step driver.

use core::any::Any;
use core::fmt;
use core::ops::Range;
use std::rc::Rc;

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{
    CarriedStepAudit, Carrier, Class, Component, ContractViolation, EndOfStep, LoadLine,
    MassStepAudit, PortLayout, PortRole, PortSpec, PortTangent, Proposal, Receipt, STREAM_CARRIED,
    STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, audit_carried_step, audit_mass_step, audit_step,
};

use crate::erased::{Node, OnLine, OwnerError, Plain, Tangent};
use crate::ledger::{BoundaryAccount, ConnectionAccount, Ledger, NodeAccount, SinkAccount, Sum};
use crate::solve::Condensed;
use crate::wiring::{
    AdmissionError, Block, Connection, Driver, HeatSink, MapMember, NodeId, PortRef,
    RouteTemperature, blocks_in_dependency_order, tolerance_group,
};

mod line;

use line::Line;

/// A prescribed port input: `(t0, dt, values)` fills the value held over the
/// step that starts at `t0` and lasts `dt`.
pub type BoundaryFn = Box<dyn Fn(f64, f64, &mut [f64])>;

/// How many units of first-order rounding a block solve accepts in a
/// residual. The estimate is `eps * sum |d residual / d input| |input|` over
/// the block's own inputs; the factor covers the few operations between an
/// input and a residual and the inputs a block does not solve for.
const ROUNDING_UNITS: f64 = 16.0;

/// A differenced entry of a column is resolved when the change it is
/// the quotient of is at least this share of its output, `2^-32`, some
/// `2^20` units in the output's last place: the output's rounding then
/// leaves the entry within about 1e-6 of itself. Numerical policy.
const RESOLVED: f64 = 1_048_576.0 * f64::EPSILON;

/// What a longer step aims an unresolved entry's change at: `2^-26` of its
/// output, `sqrt(eps)`, some `2^26` units in its last place, `2^6` times
/// what resolves it. Numerical policy.
const AIMED: f64 = 67_108_864.0 * f64::EPSILON;

/// A longer step that changes an output by more than this share of it,
/// `2^-13`, `2^13` times what the step aimed at, has gone past where the
/// output answers the input in proportion: across a jump. The entry keeps
/// what it had.
const LOCAL: f64 = 8_192.0 * AIMED;

/// The longest step a differenced column is taken at, as a share
/// of the input's scale, or of one unit where that is less: `2^-10`.
const WIDEST: f64 = 1.0 / 1_024.0;

/// How far an entry taken at a longer step may lie from the same entry
/// differenced at that step on the other side of the iterate, as a share of
/// it, and be kept: `2^-10`. Further apart, the output does not answer in
/// proportion over the two steps (it has a kink there), and the entry is
/// left as the first step had it.
const AGREED: f64 = 1.0 / 1_024.0;

/// How many times a differenced column is differenced again at a longer
/// step, at most, besides the step on the other side.
const LONGER_STEPS: usize = 4;

/// Whether a differenced entry whose output `current` changed by `change`
/// (a magnitude) is resolved: by at least [`RESOLVED`] of the output, or
/// the output is zero.
fn resolved(change: f64, current: f64) -> bool {
    current == 0.0 || change >= RESOLVED * current.abs()
}

/// How many times longer a step has to be to aim an unresolved change at
/// [`AIMED`] of its output; a change that rounded to nothing is taken
/// `2^26` times longer.
fn aimed(change: f64, current: f64) -> f64 {
    if change > 0.0 {
        AIMED * current.abs() / change
    } else {
        1.0 / AIMED
    }
}

/// Where an entry of a differenced column stands while the column is
/// resolved (`System::resolve_column`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Entry {
    /// Resolved by the first step: the entry it was.
    First,
    /// Not resolved yet: the difference at the latest step.
    Open,
    /// Resolved at a longer step, kept if the other side agrees.
    Taken,
    /// Left unresolved.
    Left,
}

/// Absolute term sum of a power pairing, before cancellation. Ordinary
/// power ports retain their original component-wise arithmetic.
fn power_pairing_scale(spec: &PortSpec, input: &[f64], output: &[f64]) -> f64 {
    if let Some(metric) = spec.field_metric() {
        let (effort, flow) = if spec.role() == PortRole::EffortIn {
            (input, output)
        } else {
            (output, input)
        };
        let mut scale = 0.0;
        for (i, &e) in effort.iter().enumerate() {
            for (j, &f) in flow.iter().enumerate() {
                scale += e.abs() * (metric.entries()[i * metric.width() + j] * f).abs();
            }
        }
        scale
    } else {
        input.iter().zip(output).map(|(e, f)| (e * f).abs()).sum()
    }
}

/// Energy allowed by one input coefficient's solve error. Both roles use
/// the same metric because it is symmetric. Apply the small allowance
/// before the metric so an unrepresentable conjugate need not overflow a
/// finite energy allowance.
fn power_input_slack_j(
    spec: &PortSpec,
    output: &[f64],
    coefficient: usize,
    dt_s: f64,
    allowance: f64,
) -> f64 {
    if allowance == 0.0 {
        return 0.0;
    }
    let scale = dt_s * allowance;
    if let Some(metric) = spec.field_metric() {
        metric.entries()[coefficient * metric.width()..(coefficient + 1) * metric.width()]
            .iter()
            .zip(output)
            .map(|(m, value)| m * (scale * value))
            .sum::<f64>()
            .abs()
    } else {
        scale * output[coefficient].abs()
    }
}

/// Limits of the step driver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveConfig {
    /// Interface residual allowed, relative to the larger of the two sides
    /// and to the largest unknown of the same kind in the block.
    ///
    /// The same tolerance is applied a second time to what the residuals
    /// do: carried once through the components that read them, by the
    /// slopes of the block's Jacobian formed in the same step, they may
    /// not move any wire's target by more than that wire's own tolerance.
    /// A temperature
    /// is thereby closed to what the heat it drives requires, not to a
    /// fraction of its absolute value. Only the part of a residual above
    /// its own rounding floor is carried.
    pub relative_tolerance: f64,
    /// Interface residual always allowed, in the unknown's own unit.
    pub absolute_tolerance: f64,
    /// Newton iterations allowed per block and trial.
    pub max_iterations: u32,
    /// Times one step may be shortened at a component's request.
    pub max_shortenings: u32,
    /// Interface energy a connection may fail to cancel, relative to the
    /// energy that crossed it.
    pub interface_energy_tolerance: f64,
}

impl Default for SolveConfig {
    fn default() -> Self {
        Self {
            relative_tolerance: 1.0e-12,
            absolute_tolerance: 0.0,
            max_iterations: 25,
            max_shortenings: 64,
            interface_energy_tolerance: 1.0e-9,
        }
    }
}

/// Where in a step a component refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// `limit_dt_s`.
    Negotiate,
    /// `propose`, or the tangent of its port map.
    Propose,
    /// Stored energy or mass measurement during the audit.
    Measure,
    /// `admit`.
    Admit,
}

/// Why a step was not taken. Nothing is published when a step fails.
#[derive(Debug)]
pub enum SystemError {
    /// A connection failed to cancel a carried stream coordinate.
    InterfaceCarried {
        /// Connection in declaration order.
        connection: usize,
        /// Zero-based carried coordinate.
        carried: usize,
        /// Uncancelled transported amount.
        defect: f64,
        /// Allowance in the coordinate's declared unit.
        tolerance: f64,
    },
    /// A component refused with its own error.
    Component {
        /// Component name.
        node: String,
        /// Where it refused.
        phase: Phase,
        /// The owner's error, unchanged.
        error: OwnerError,
    },
    /// A component's step broke the contract.
    Contract {
        /// Component name.
        node: String,
        /// The law it broke.
        violation: ContractViolation,
    },
    /// A component answered a negotiation outside the contract.
    Negotiation {
        /// Component name.
        node: String,
        /// Step offered to it [s].
        proposed_s: f64,
        /// Step it returned [s].
        returned_s: f64,
    },
    /// A block's interface equations did not converge.
    NotConverged {
        /// Components of the block.
        nodes: Vec<String>,
        /// Iterations spent.
        iterations: u32,
        /// As `node.port`: the input whose residual is furthest beyond its
        /// own tolerance; or, when every residual passes its own test, the
        /// open input that moves another wire's target the most.
        worst: String,
        /// Its residual.
        residual: f64,
    },
    /// A block's interface Jacobian is singular.
    Singular {
        /// Components of the block.
        nodes: Vec<String>,
    },
    /// A component inside a loop offers no tangent and is not smooth.
    Unsupported {
        /// Component name.
        node: String,
        /// What the driver would need.
        needs: &'static str,
    },
    /// The load line of a component could not be used.
    LoadLine {
        /// Component name.
        node: String,
        /// Why.
        reason: &'static str,
    },
    /// A temperature inside a loop has no positive starting value.
    GuessRequired {
        /// The input, as `node.port`.
        port: String,
    },
    /// A connection failed to cancel interface energy within tolerance.
    Interface {
        /// Connection index, in declaration order.
        connection: usize,
        /// Energy that failed to cancel [J].
        defect_j: f64,
        /// Tolerance applied [J].
        tolerance_j: f64,
    },
    /// A connection failed to cancel stream mass within its solve allowance.
    InterfaceMass {
        /// Connection index, in declaration order.
        connection: usize,
        /// Mass that failed to cancel [kg].
        defect_kg: f64,
        /// Tolerance applied [kg].
        tolerance_kg: f64,
    },
    /// The step would destroy entropy.
    SecondLaw {
        /// Entropy generated in the step [J/K].
        generated_j_per_k: f64,
    },
    /// The step was shortened more often than the configuration allows.
    StepUnderflow {
        /// Last step tried [s].
        dt_s: f64,
    },
    /// The requested step is not positive and finite.
    Step {
        /// The rejected step [s].
        dt_s: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for SystemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InterfaceCarried { connection, carried, defect, tolerance } => write!(f,
                "connection {connection} left carried coordinate {carried} uncancelled by {defect} (tolerance {tolerance})"),
            Self::Component { node, phase, error } => {
                write!(f, "component `{node}` refused in {phase:?}: {error:?}")
            }
            Self::Contract { node, violation } => write!(f, "component `{node}`: {violation}"),
            Self::Negotiation {
                node,
                proposed_s,
                returned_s,
            } => write!(
                f,
                "component `{node}` answered {returned_s} s to a proposed step of {proposed_s} s"
            ),
            Self::NotConverged {
                nodes,
                iterations,
                worst,
                residual,
            } => write!(
                f,
                "block {nodes:?} did not converge in {iterations} iterations; \
                 worst input `{worst}` has residual {residual}"
            ),
            Self::Singular { nodes } => {
                write!(f, "block {nodes:?} has a singular interface Jacobian")
            }
            Self::Unsupported { node, needs } => {
                write!(f, "component `{node}` is inside a loop and needs {needs}")
            }
            Self::LoadLine { node, reason } => {
                write!(f, "component `{node}` cannot be put on its load line: {reason}")
            }
            Self::GuessRequired { port } => write!(
                f,
                "temperature input `{port}` is inside a loop and needs a positive starting value"
            ),
            Self::Interface {
                connection,
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "connection {connection} left {defect_j} J uncancelled (tolerance {tolerance_j} J)"
            ),
            Self::InterfaceMass {
                connection,
                defect_kg,
                tolerance_kg,
            } => write!(
                f,
                "connection {connection} left {defect_kg} kg uncancelled (tolerance {tolerance_kg} kg)"
            ),
            Self::SecondLaw { generated_j_per_k } => {
                write!(
                    f,
                    "the step would generate {generated_j_per_k} J/K of entropy"
                )
            }
            Self::StepUnderflow { dt_s } => {
                write!(f, "the step was shortened too often (last {dt_s} s)")
            }
            Self::Step { dt_s } => write!(f, "step {dt_s} s is not positive and finite"),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for SystemError {}

/// Why an accepted-state checkpoint could not be captured or restored.
#[derive(Debug)]
pub enum CheckpointError {
    /// The checkpoint belongs to another sealed system, even if its graph
    /// has the same component names and types. Nothing was restored.
    ForeignSystem,
    /// An erased state slot does not match its owner's concrete state type.
    /// This is a kernel invariant violation, not a component refusal.
    State {
        /// Component name.
        node: String,
        /// The type-erasure error, recoverable with `downcast_ref`.
        source: OwnerError,
    },
}

impl fmt::Display for CheckpointError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignSystem => write!(f, "checkpoint belongs to another sealed system"),
            Self::State { node, source } => {
                write!(
                    f,
                    "component `{node}` cannot clone its accepted state: {source:?}"
                )
            }
        }
    }
}

impl std::error::Error for CheckpointError {}

/// A reusable in-memory checkpoint of one sealed [`System`].
///
/// It retains accepted component states, interface values and warm starts,
/// the compensated clock, step count, complete ledger and solve configuration.
/// It can be inspected after its system is dropped, but can be restored only
/// into the same system that created it. It contains no candidate trial state.
///
/// This is not an encoded checkpoint for process restart or graph cloning.
/// Boundaries and components remain those of the sealed system. Any external
/// data captured by a [`BoundaryFn`] is outside this checkpoint; bitwise replay
/// requires the same boundary responses and pure component proposals.
pub struct Checkpoint {
    identity: Rc<()>,
    accepted: Vec<Box<dyn Any>>,
    inputs: Vec<f64>,
    outputs: Vec<f64>,
    time: Sum,
    steps: u64,
    ledger: Ledger,
    config: SolveConfig,
}

impl fmt::Debug for Checkpoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Checkpoint")
            .field("components", &self.accepted.len())
            .field("time_s", &self.time.value())
            .field("steps", &self.steps)
            .field("config", &self.config)
            .finish_non_exhaustive()
    }
}

impl Checkpoint {
    /// Simulated time captured by the checkpoint [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.time.value()
    }

    /// Accepted steps captured by the checkpoint.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// The saved accepted state of `node`, if it has type `S`.
    /// Use the node identifier from the originating system.
    #[must_use]
    pub fn state<S: 'static>(&self, node: NodeId) -> Option<&S> {
        self.accepted.get(node.0)?.downcast_ref::<S>()
    }

    /// The complete window ledger, including each sum's compensation.
    #[must_use]
    pub fn ledger(&self) -> &Ledger {
        &self.ledger
    }

    /// The solve configuration restored with this checkpoint.
    #[must_use]
    pub fn config(&self) -> SolveConfig {
        self.config
    }
}

/// What one accepted step did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepReport {
    /// Step taken [s].
    pub dt_s: f64,
    /// Time after the step [s].
    pub time_s: f64,
    /// Newton iterations summed over blocks.
    pub iterations: u32,
    /// Times the step was shortened before it was taken.
    pub shortenings: u32,
    /// Entries of differenced columns, a load line's or a block solve's,
    /// left unresolved against the rounding of their outputs, over the
    /// trials of the step (the contract, "Differenced columns"). The step is
    /// taken and its balance audited all the same: the count says that a
    /// line component's relation was solved, or an update made, on a matrix
    /// some of whose entries are a count of rounding units, the first step's
    /// difference, or nothing. An output that does not answer an input at
    /// all (a port of a two-port neighbour that the force at its other port
    /// does not move) is counted at every step, and its entry, nothing, is
    /// right. Above that, a caller that needs the line component's answer
    /// to more than the step's balance (whether a pair sticks or slips, and
    /// when) declares the neighbour's tangent
    /// ([`SystemBuilder::add_with_tangent`]), which the line is handed as
    /// declared.
    pub unresolved: u32,
}

/// Collects components and connections, then seals them into a [`System`].
#[derive(Default)]
pub struct SystemBuilder {
    names: Vec<String>,
    nodes: Vec<Box<dyn Node>>,
    layouts: Vec<PortLayout>,
    accepted: Vec<Box<dyn Any>>,
    candidate: Vec<Box<dyn Any>>,
    used: Vec<Vec<bool>>,
    /// Ports that receive heat from a heat route; another route may add to
    /// them.
    route_receiver: Vec<Vec<bool>>,
    connections: Vec<Connection>,
    boundaries: Vec<(PortRef, BoundaryFn)>,
    sinks: Vec<HeatSink>,
}

impl SystemBuilder {
    /// An empty description.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn insert(
        &mut self,
        name: &str,
        node: Box<dyn Node>,
        accepted: Box<dyn Any>,
        candidate: Box<dyn Any>,
    ) -> Result<NodeId, AdmissionError> {
        if self.names.iter().any(|known| known == name) {
            return Err(AdmissionError::DuplicateNode {
                name: name.to_owned(),
            });
        }
        let layout = PortLayout::new(node.ports()).map_err(|violation| AdmissionError::Layout {
            node: name.to_owned(),
            violation,
        })?;
        self.used.push(vec![false; node.ports().len()]);
        self.route_receiver.push(vec![false; node.ports().len()]);
        self.names.push(name.to_owned());
        self.layouts.push(layout);
        self.nodes.push(node);
        self.accepted.push(accepted);
        self.candidate.push(candidate);
        Ok(NodeId(self.names.len() - 1))
    }

    /// Add a component with its initial state.
    ///
    /// # Errors
    /// A repeated name or an invalid port list.
    pub fn add<C>(
        &mut self,
        name: &str,
        component: C,
        state: C::State,
    ) -> Result<NodeId, AdmissionError>
    where
        C: Component + 'static,
        C::State: 'static,
    {
        let candidate = Box::new(state.clone());
        self.insert(name, Box::new(Plain(component)), Box::new(state), candidate)
    }

    /// Add a component that offers the tangent of its port map.
    ///
    /// # Errors
    /// A repeated name or an invalid port list.
    pub fn add_with_tangent<C>(
        &mut self,
        name: &str,
        component: C,
        state: C::State,
    ) -> Result<NodeId, AdmissionError>
    where
        C: PortTangent + 'static,
        C::State: 'static,
    {
        let candidate = Box::new(state.clone());
        self.insert(
            name,
            Box::new(Tangent(component)),
            Box::new(state),
            candidate,
        )
    }

    /// Add a component whose relation at its power ports is not a function
    /// of the flow, to be stepped on the load line its neighbours present
    /// ([`LoadLine`]): dry friction, a contact.
    ///
    /// Inside a loop the solve finds the line once per step and the
    /// component settles on it; one such component per loop. Where every
    /// flow of its line is prescribed or comes from an earlier block, there
    /// is no line to find and the component is stepped by its mandatory
    /// tier.
    ///
    /// # Errors
    /// A repeated name, an invalid port list, or a line that is not over
    /// power ports that receive a flow, in increasing order.
    pub fn add_on_line<C>(
        &mut self,
        name: &str,
        component: C,
        state: C::State,
    ) -> Result<NodeId, AdmissionError>
    where
        C: LoadLine + 'static,
        C::State: 'static,
    {
        let ports = component.ports();
        let line = component.line_ports();
        let increasing = line.windows(2).all(|pair| pair[0] < pair[1]);
        let receive_flow = line.iter().all(|&port| {
            ports.get(port).is_some_and(|spec| {
                matches!(spec.carrier(), Carrier::Power(_)) && spec.role() == PortRole::FlowIn
            })
        });
        if line.is_empty() || !increasing || !receive_flow {
            return Err(AdmissionError::LinePorts {
                node: name.to_owned(),
            });
        }
        let candidate = Box::new(state.clone());
        self.insert(name, Box::new(OnLine(component)), Box::new(state), candidate)
    }

    /// The port of `node` with identifier `id`.
    ///
    /// # Errors
    /// [`AdmissionError::UnknownPort`].
    pub fn port(&self, node: NodeId, id: &str) -> Result<PortRef, AdmissionError> {
        self.nodes
            .get(node.0)
            .and_then(|entry| {
                entry
                    .ports()
                    .iter()
                    .position(|port| port.id().as_str() == id)
            })
            .map(|port| PortRef { node: node.0, port })
            .ok_or_else(|| AdmissionError::UnknownPort {
                node: self.names.get(node.0).cloned().unwrap_or_default(),
                port: id.to_owned(),
            })
    }

    fn spec(&self, port: PortRef) -> &PortSpec {
        &self.nodes[port.node].ports()[port.port]
    }

    fn label(&self, port: PortRef) -> String {
        format!(
            "{}.{}",
            self.names[port.node],
            self.spec(port).id().as_str()
        )
    }

    /// Refuse a port that is already connected or named twice in one call.
    ///
    /// Every connection method checks all of its ports with this before it
    /// changes anything, so a refused call leaves the builder as it was.
    fn check_free(&self, ports: &[PortRef]) -> Result<(), AdmissionError> {
        for (index, &port) in ports.iter().enumerate() {
            if self.used[port.node][port.port] || ports[..index].contains(&port) {
                return Err(AdmissionError::AlreadyConnected {
                    port: self.label(port),
                });
            }
        }
        Ok(())
    }

    fn take(&mut self, ports: &[PortRef]) {
        for &port in ports {
            self.used[port.node][port.port] = true;
        }
    }

    fn need_role(
        &self,
        port: PortRef,
        role: PortRole,
        needs: &'static str,
    ) -> Result<(), AdmissionError> {
        if self.spec(port).role() == role {
            Ok(())
        } else {
            Err(AdmissionError::Role {
                port: self.label(port),
                needs,
            })
        }
    }

    fn need_match(&self, a: PortRef, b: PortRef) -> Result<(), AdmissionError> {
        let (left, right) = (self.spec(a), self.spec(b));
        if left.carrier() == right.carrier() && left.width() == right.width() {
            if left.carrier() == Carrier::Stream && left.stream_basis() != right.stream_basis() {
                return Err(AdmissionError::StreamBasis {
                    a: self.label(a),
                    b: self.label(b),
                });
            }
            self.need_pairing_match(a, b)
        } else {
            Err(AdmissionError::Mismatch {
                a: self.label(a),
                b: self.label(b),
            })
        }
    }

    /// Direct copies and scalar ratios preserve power only in the same
    /// pairing. An explicit identity also matches an ordinary power port.
    fn need_pairing_match(&self, a: PortRef, b: PortRef) -> Result<(), AdmissionError> {
        let equal = match (self.spec(a).field_metric(), self.spec(b).field_metric()) {
            (None, None) => true,
            (Some(left), Some(right)) => left == right,
            (Some(metric), None) | (None, Some(metric)) => {
                metric.entries().iter().enumerate().all(|(index, &value)| {
                    value
                        == if index / metric.width() == index % metric.width() {
                            1.0
                        } else {
                            0.0
                        }
                })
            }
        };
        if equal {
            Ok(())
        } else {
            Err(AdmissionError::Mismatch {
                a: self.label(a),
                b: self.label(b),
            })
        }
    }

    /// Admit a junction. A plain junction joins ports that carry the same
    /// thing; a `geared` one joins power ports of any domain, each member
    /// through its own nonzero ratio.
    fn junction(
        &mut self,
        setter: PortRef,
        members: &[(PortRef, f64)],
        setter_role: PortRole,
        member_role: PortRole,
        geared: bool,
    ) -> Result<(), AdmissionError> {
        if members.is_empty() {
            return Err(AdmissionError::Empty);
        }
        if geared {
            self.need_power(setter, setter_role)?;
            for &(member, _) in members {
                self.need_power(member, member_role)?;
                if self.spec(setter).width() != self.spec(member).width() {
                    return Err(AdmissionError::Mismatch {
                        a: self.label(setter),
                        b: self.label(member),
                    });
                }
                self.need_pairing_match(setter, member)?;
            }
            if let Some(&(_, ratio)) = members
                .iter()
                .find(|(_, ratio)| !ratio.is_finite() || *ratio == 0.0)
            {
                return Err(AdmissionError::TransformerRatio { ratio });
            }
        } else {
            let (setter_needs, member_needs) = if setter_role == PortRole::FlowIn {
                ("a port that receives flow", "a port that receives effort")
            } else {
                ("a port that receives effort", "a port that receives flow")
            };
            self.need_role(setter, setter_role, setter_needs)?;
            for &(member, _) in members {
                self.need_role(member, member_role, member_needs)?;
                self.need_match(setter, member)?;
            }
        }
        // Heat rate is power itself, so handing the same heat rate to several
        // members would multiply the energy. Heat ports share a temperature
        // (common effort) or go through a heat route.
        // A stream's rates are what crosses, in the same way: several
        // storages cannot each be handed the whole of one flow, and a
        // transport element has no use for the sum of their pressures. A
        // storage is the node that streams meet at (common effort).
        if setter_role == PortRole::EffortIn && members.len() > 1 {
            match self.spec(setter).carrier() {
                Carrier::Heat => {
                    return Err(AdmissionError::HeatCommonFlow {
                        port: self.label(setter),
                    });
                }
                Carrier::Stream => {
                    return Err(AdmissionError::StreamCommonFlow {
                        port: self.label(setter),
                    });
                }
                Carrier::Power(_) | Carrier::Signal => {}
            }
        }
        let mut all = Vec::with_capacity(members.len() + 1);
        all.push(setter);
        all.extend(members.iter().map(|(port, _)| *port));
        self.check_free(&all)?;
        self.take(&all);
        Ok(())
    }

    /// Join ports at a common effort: `setter` returns the effort and receives
    /// minus the sum of the members' flows.
    ///
    /// # Errors
    /// Wrong roles, mismatched ports, or a port already connected.
    pub fn common_effort(
        &mut self,
        setter: PortRef,
        members: &[PortRef],
    ) -> Result<(), AdmissionError> {
        let members: Vec<_> = members.iter().map(|&port| (port, 1.0)).collect();
        self.junction(
            setter,
            &members,
            PortRole::FlowIn,
            PortRole::EffortIn,
            false,
        )?;
        self.connections
            .push(Connection::CommonEffort { setter, members });
        Ok(())
    }

    /// Join ports at a common flow: `setter` returns the flow, every member
    /// receives minus that flow, and the setter receives the sum of the
    /// members' efforts.
    ///
    /// # Errors
    /// Wrong roles, mismatched ports, or a port already connected.
    pub fn common_flow(
        &mut self,
        setter: PortRef,
        members: &[PortRef],
    ) -> Result<(), AdmissionError> {
        let members: Vec<_> = members.iter().map(|&port| (port, 1.0)).collect();
        self.junction(
            setter,
            &members,
            PortRole::EffortIn,
            PortRole::FlowIn,
            false,
        )?;
        self.connections
            .push(Connection::CommonFlow { setter, members });
        Ok(())
    }

    /// Join power ports at a common flow, each member through its own ideal
    /// transformer.
    ///
    /// `setter` receives effort and returns the flow `f`. Member `i` receives
    /// flow `-ratio_i f` and returns an effort `e_i`; the setter receives
    /// `sum(ratio_i e_i)`. The power entering all owners cancels for any
    /// ratios. One member is [`Self::transformer`]; every ratio equal to one
    /// is [`Self::common_flow`] between power ports.
    ///
    /// The members may be in different power domains: each ratio supplies the
    /// conversion to the setter's in the ports' SI units, and its sign the
    /// orientation. A wall of area `A` between two gases, pushed by a force,
    /// is one junction: the wall's inertia sets the velocity, one gas has
    /// ratio `A`, the other `-A`, and the force has ratio one.
    ///
    /// # Errors
    /// Wrong roles, a heat or signal port, unequal widths, a zero or
    /// non-finite ratio, no member, or a port already connected. Refusal
    /// leaves the builder unchanged.
    pub fn common_flow_ratios(
        &mut self,
        setter: PortRef,
        members: &[(PortRef, f64)],
    ) -> Result<(), AdmissionError> {
        self.junction(setter, members, PortRole::EffortIn, PortRole::FlowIn, true)?;
        self.connections.push(Connection::CommonFlow {
            setter,
            members: members.to_vec(),
        });
        Ok(())
    }

    /// Join power ports at a common effort, each member through its own ideal
    /// transformer: the dual of [`Self::common_flow_ratios`].
    ///
    /// `setter` receives flow and returns the effort `e`. Member `i` receives
    /// effort `ratio_i e` and returns a flow `f_i`; the setter receives
    /// `-sum(ratio_i f_i)`. The power entering all owners cancels for any
    /// ratios.
    ///
    /// # Errors
    /// As for [`Self::common_flow_ratios`].
    pub fn common_effort_ratios(
        &mut self,
        setter: PortRef,
        members: &[(PortRef, f64)],
    ) -> Result<(), AdmissionError> {
        self.junction(setter, members, PortRole::FlowIn, PortRole::EffortIn, true)?;
        self.connections.push(Connection::CommonEffort {
            setter,
            members: members.to_vec(),
        });
        Ok(())
    }

    /// Join two ports of the same kind and opposite roles.
    ///
    /// # Errors
    /// Roles that are not one effort-receiving and one flow-receiving port,
    /// mismatched ports, or a port already connected.
    pub fn bond(&mut self, a: PortRef, b: PortRef) -> Result<(), AdmissionError> {
        if self.spec(a).role() == PortRole::FlowIn {
            self.common_effort(a, &[b])
        } else {
            self.common_effort(b, &[a])
        }
    }

    /// Require a power port with `role`, which is effort- or flow-receiving.
    fn need_power(&self, port: PortRef, role: PortRole) -> Result<(), AdmissionError> {
        let needs = if role == PortRole::EffortIn {
            "a power port that receives effort"
        } else {
            "a power port that receives flow"
        };
        if !matches!(self.spec(port).carrier(), Carrier::Power(_)) {
            return Err(AdmissionError::Role {
                port: self.label(port),
                needs,
            });
        }
        self.need_role(port, role, needs)
    }

    fn need_power_map_roles(&self, a: PortRef, b: PortRef) -> Result<(), AdmissionError> {
        self.need_power(a, PortRole::EffortIn)?;
        self.need_power(b, PortRole::FlowIn)
    }

    /// Connect power ports through an ideal transformer of signed `ratio`.
    ///
    /// `a` receives effort and returns flow; `b` receives flow and returns
    /// effort. The generated equations are `f_b = -ratio f_a` and
    /// `e_a = ratio e_b`, so the power entering the two owners cancels.
    /// The ratio acts on each coefficient when the ports have width > 1.
    ///
    /// Unlike a bond, this explicit transformation may join different power
    /// domains: a lever, gear, piston area or electrical turns ratio supplies
    /// the conversion in the ports' SI units. A negative ratio reverses the
    /// orientation. Heat-rate and signal ports are not power-conjugate pairs
    /// and are refused. Two integral-causality inertias cannot be joined by
    /// this connection: one side must receive flow.
    ///
    /// # Errors
    /// Wrong roles or carriers, unequal widths, a zero or non-finite ratio,
    /// or a port already connected. Refusal leaves the builder unchanged.
    pub fn transformer(
        &mut self,
        a: PortRef,
        b: PortRef,
        ratio: f64,
    ) -> Result<(), AdmissionError> {
        self.need_power_map_roles(a, b)?;
        if self.spec(a).width() != self.spec(b).width() {
            return Err(AdmissionError::Mismatch {
                a: self.label(a),
                b: self.label(b),
            });
        }
        self.need_pairing_match(a, b)?;
        if !ratio.is_finite() || ratio == 0.0 {
            return Err(AdmissionError::TransformerRatio { ratio });
        }
        self.check_free(&[a, b])?;
        self.take(&[a, b]);
        self.connections
            .push(Connection::Transformer { a, b, ratio });
        Ok(())
    }

    /// Connect power ports of possibly different widths by an adjoint pair.
    ///
    /// `a` receives effort; `b` receives flow. `map` is the row-major matrix
    /// `P` with `b.width()` rows and `a.width()` columns. The kernel generates
    /// `f_b = -P f_a` and `e_a = M_a^-1 P^T M_b e_b`, so the declared
    /// power pairings cancel, including for rectangular or rank-deficient
    /// maps. Ordinary ports use the identity metric. The effort map is
    /// computed at admission using the metrics' cached Cholesky factors;
    /// a step performs no metric solve. Different power domains are allowed;
    /// coefficients must supply the conversion in the ports' SI units.
    /// This does not infer a trace metric or a conservative heat map.
    ///
    /// This is [`Self::power_map_members`] with one member.
    ///
    /// # Errors
    /// Wrong roles or carriers, a matrix of the wrong size, a non-finite
    /// coefficient, an unrepresentable metric adjoint, or a port already
    /// connected. Refusal changes nothing.
    pub fn power_map(
        &mut self,
        a: PortRef,
        b: PortRef,
        map: Vec<f64>,
    ) -> Result<(), AdmissionError> {
        self.power_map_members(vec![(a, map)], b)
    }

    /// Connect several effort-receiving power ports to one flow-receiving
    /// port by one adjoint pair.
    ///
    /// Member `k` is `(a_k, P_k)`, with `P_k` row-major of `b.width()` rows
    /// and `a_k.width()` columns. The kernel generates `f_b = -sum(P_k f_a_k)`
    /// and `e_a_k = M_a_k^-1 P_k^T M_b e_b`, so the declared power pairings
    /// cancel. It is the map `[P_1 ... P_n]` of the members' coefficients
    /// taken together: a rigid lever from the force and torque ports of two
    /// bodies onto one contact, for instance. With one member it is
    /// [`Self::power_map`]. [`Self::common_effort_ratios`] is its scalar
    /// case: a member with ratio `r` is `P = r I`.
    ///
    /// The map is constant from admission. A lever whose arms turn with its
    /// bodies keeps the ratios of the pose it was admitted at.
    ///
    /// # Errors
    /// No member ([`AdmissionError::Empty`]), and, member by member, what
    /// [`Self::power_map`] refuses; the index of a non-finite coefficient
    /// counts through the members' maps in order. A port named twice is
    /// already connected. Refusal changes nothing.
    pub fn power_map_members(
        &mut self,
        members: Vec<(PortRef, Vec<f64>)>,
        b: PortRef,
    ) -> Result<(), AdmissionError> {
        if members.is_empty() {
            return Err(AdmissionError::Empty);
        }
        let mut before = 0;
        for (a, map) in &members {
            self.need_power_map_roles(*a, b)?;
            let (rows, columns) = (self.spec(b).width(), self.spec(*a).width());
            if rows.checked_mul(columns) != Some(map.len()) {
                return Err(AdmissionError::MapShape {
                    rows,
                    columns,
                    entries: map.len(),
                });
            }
            if let Some((index, &value)) = map.iter().enumerate().find(|(_, v)| !v.is_finite()) {
                return Err(AdmissionError::MapCoefficient {
                    index: before + index,
                    value,
                });
            }
            before += map.len();
        }
        let mut ports: Vec<PortRef> = members.iter().map(|(a, _)| *a).collect();
        ports.push(b);
        self.check_free(&ports)?;
        let mut mapped = Vec::with_capacity(members.len());
        for (a, map) in members {
            let effort_map = self.field_effort_map(a, b, &map)?;
            mapped.push(MapMember { a, map, effort_map });
        }
        self.take(&ports);
        self.connections
            .push(Connection::PowerMap { members: mapped, b });
        Ok(())
    }

    /// Construct the metric adjoint once; ordinary ports keep their exact
    /// existing transpose path and do not allocate another matrix.
    fn field_effort_map(
        &self,
        a: PortRef,
        b: PortRef,
        map: &[f64],
    ) -> Result<Option<Vec<f64>>, AdmissionError> {
        let (left, right) = (self.spec(a), self.spec(b));
        if left.field_metric().is_none() && right.field_metric().is_none() {
            return Ok(None);
        }
        let (columns, rows) = (left.width(), right.width());
        let mut result = vec![0.0; map.len()];
        let mut rhs = vec![0.0; columns];
        let mut solved = vec![0.0; columns];
        let refusal = |violation| AdmissionError::MapAdjoint {
            a: self.label(a),
            b: self.label(b),
            violation,
        };
        for j in 0..rows {
            for i in 0..columns {
                rhs[i] = if let Some(metric) = right.field_metric() {
                    (0..rows)
                        .map(|k| map[k * columns + i] * metric.entries()[k * rows + j])
                        .sum()
                } else {
                    map[j * columns + i]
                };
            }
            if rhs.iter().any(|value| !value.is_finite()) {
                return Err(refusal(ContractViolation::NonFinite {
                    quantity: "field effort adjoint",
                }));
            }
            if let Some(metric) = left.field_metric() {
                metric.solve(&rhs, &mut solved).map_err(&refusal)?;
            } else {
                solved.copy_from_slice(&rhs);
            }
            for i in 0..columns {
                result[i * rows + j] = solved[i];
            }
        }
        Ok(Some(result))
    }

    /// Deliver the heat leaving `source` to `sinks` by fraction. Fractions
    /// below one in total need `remainder`, a named sink outside the system.
    ///
    /// A receiver may be named by several routes: it receives the sum of
    /// their shares, and those routes are audited together. A source belongs
    /// to one route, and a receiver of a route takes no other connection.
    ///
    /// # Errors
    /// Wrong roles or carriers, bad fractions, a missing sink, a bad
    /// temperature, a port named twice, or a port already connected to
    /// anything but another route's delivery.
    pub fn heat_route(
        &mut self,
        source: PortRef,
        sinks: &[(PortRef, f64)],
        remainder: Option<HeatSink>,
        temperature: RouteTemperature,
    ) -> Result<(), AdmissionError> {
        let heat_source = "a heat port that receives temperature";
        let heat_sink = "a heat port that receives heat rate";
        if self.spec(source).carrier() != Carrier::Heat {
            return Err(AdmissionError::Role {
                port: self.label(source),
                needs: heat_source,
            });
        }
        self.need_role(source, PortRole::EffortIn, heat_source)?;
        let mut sum = 0.0;
        for &(sink, fraction) in sinks {
            if self.spec(sink).carrier() != Carrier::Heat {
                return Err(AdmissionError::Role {
                    port: self.label(sink),
                    needs: heat_sink,
                });
            }
            self.need_role(sink, PortRole::FlowIn, heat_sink)?;
            if !(fraction > 0.0 && fraction <= 1.0) {
                return Err(AdmissionError::Fractions { sum: fraction });
            }
            sum += fraction;
        }
        let slack = 8.0 * f64::EPSILON * (sinks.len() as f64 + 1.0);
        if !(sum <= 1.0 + slack) {
            return Err(AdmissionError::Fractions { sum });
        }
        let unallocated = (1.0 - sum).max(0.0);
        if let RouteTemperature::Fixed(kelvin) = temperature
            && !(kelvin.is_finite() && kelvin > 0.0)
        {
            return Err(AdmissionError::Temperature { kelvin });
        }
        match &remainder {
            Some(sink) if !(sink.temperature_k.is_finite() && sink.temperature_k > 0.0) => {
                return Err(AdmissionError::Temperature {
                    kelvin: sink.temperature_k,
                });
            }
            None if unallocated > slack => {
                return Err(AdmissionError::SinkRequired { unallocated });
            }
            None if sinks.is_empty() => return Err(AdmissionError::Empty),
            _ => {}
        }
        // A sink named before is the same place: one ledger row, and so one
        // temperature.
        let known = remainder
            .as_ref()
            .and_then(|sink| self.sinks.iter().position(|known| known.name == sink.name));
        if let (Some(position), Some(sink)) = (known, &remainder)
            && self.sinks[position].temperature_k.to_bits() != sink.temperature_k.to_bits()
        {
            return Err(AdmissionError::SinkTemperature {
                name: sink.name.as_str().to_owned(),
                declared_k: self.sinks[position].temperature_k,
                repeated_k: sink.temperature_k,
            });
        }
        let mut all = Vec::with_capacity(sinks.len() + 1);
        all.push(source);
        all.extend(sinks.iter().map(|(sink, _)| *sink));
        for (index, &port) in all.iter().enumerate() {
            let shared = index > 0 && self.route_receiver[port.node][port.port];
            if (self.used[port.node][port.port] && !shared) || all[..index].contains(&port) {
                return Err(AdmissionError::AlreadyConnected {
                    port: self.label(port),
                });
            }
        }
        // Every check has passed: only now does the builder change.
        self.take(&all);
        for &(sink, _) in sinks {
            self.route_receiver[sink.node][sink.port] = true;
        }
        let remainder = remainder.map(|sink| {
            let temperature_k = sink.temperature_k;
            let position = known.unwrap_or_else(|| {
                self.sinks.push(sink);
                self.sinks.len() - 1
            });
            (position, unallocated, temperature_k)
        });
        self.connections.push(Connection::HeatRoute {
            source,
            sinks: sinks.to_vec(),
            remainder,
            temperature,
        });
        Ok(())
    }

    /// Copy a signal output to a signal input. One output may feed several
    /// inputs.
    ///
    /// # Errors
    /// Wrong roles, mismatched widths, or an input already driven.
    pub fn signal(&mut self, from: PortRef, to: PortRef) -> Result<(), AdmissionError> {
        self.need_role(from, PortRole::SignalOut, "a signal output")?;
        self.need_role(to, PortRole::SignalIn, "a signal input")?;
        self.need_match(from, to)?;
        self.check_free(&[to])?;
        self.take(&[to]);
        self.connections.push(Connection::Signal { from, to });
        Ok(())
    }

    /// Prescribe the input side of `port` from outside the system.
    ///
    /// # Errors
    /// A port with no input side, or a port already connected.
    pub fn boundary(&mut self, port: PortRef, value: BoundaryFn) -> Result<(), AdmissionError> {
        if self.spec(port).input_width() == 0 {
            return Err(AdmissionError::Role {
                port: self.label(port),
                needs: "a port with an input side",
            });
        }
        self.check_free(&[port])?;
        self.take(&[port]);
        self.boundaries.push((port, value));
        Ok(())
    }

    /// Check the description and build the system.
    ///
    /// # Errors
    /// [`AdmissionError::Unconnected`] for a port that exchanges power or
    /// needs an input and is connected to nothing.
    #[allow(clippy::too_many_lines)]
    pub fn seal(self) -> Result<System, AdmissionError> {
        let count = self.nodes.len();
        let mut input_start = Vec::with_capacity(count + 1);
        let mut output_start = Vec::with_capacity(count + 1);
        let mut port_start = Vec::with_capacity(count + 1);
        let (mut inputs, mut outputs, mut ports) = (0_usize, 0_usize, 0_usize);
        for (node, layout) in self.nodes.iter().zip(&self.layouts) {
            input_start.push(inputs);
            output_start.push(outputs);
            port_start.push(ports);
            inputs += layout.input_len();
            outputs += layout.output_len();
            ports += node.ports().len();
        }
        input_start.push(inputs);
        output_start.push(outputs);
        port_start.push(ports);

        for (node, used) in self.used.iter().enumerate() {
            for (port, &taken) in used.iter().enumerate() {
                let spec = &self.nodes[node].ports()[port];
                // The end inputs of a load line are wired below, by the
                // kernel.
                let wired_below = self.nodes[node]
                    .line_ends()
                    .iter()
                    .flatten()
                    .any(|&end| end == port);
                if !taken && !wired_below && spec.role() != PortRole::SignalOut {
                    return Err(AdmissionError::Unconnected {
                        port: self.label(PortRef { node, port }),
                    });
                }
            }
        }

        let input_of = |port: PortRef| -> Range<usize> {
            let range = self.layouts[port.node].input_range(port.port);
            input_start[port.node] + range.start..input_start[port.node] + range.end
        };
        let output_of = |port: PortRef| -> Range<usize> {
            let range = self.layouts[port.node].output_range(port.port);
            output_start[port.node] + range.start..output_start[port.node] + range.end
        };

        let mut drivers = vec![Driver::Unset; inputs];
        for connection in &self.connections {
            match connection {
                Connection::CommonEffort { setter, members }
                | Connection::CommonFlow { setter, members } => {
                    let common_effort = matches!(connection, Connection::CommonEffort { .. });
                    // Common effort: members copy the setter, the setter gets
                    // minus the members' sum. Common flow: members get minus
                    // the setter, the setter gets the members' sum. A member's
                    // ratio multiplies both of its directions, so the powers
                    // cancel; a ratio of one leaves the coefficients exact.
                    let (to_member, to_setter) = if common_effort {
                        (1.0, -1.0)
                    } else {
                        (-1.0, 1.0)
                    };
                    let setter_out = output_of(*setter);
                    for (offset, index) in input_of(*setter).enumerate() {
                        let terms = members
                            .iter()
                            .map(|(member, ratio)| {
                                (output_of(*member).start + offset, to_setter * ratio)
                            })
                            .collect();
                        drivers[index] = Driver::Linear(terms);
                    }
                    for (member, ratio) in members {
                        for (offset, index) in input_of(*member).enumerate() {
                            drivers[index] = Driver::Linear(vec![(
                                setter_out.start + offset,
                                to_member * ratio,
                            )]);
                        }
                    }
                }
                Connection::Transformer { a, b, ratio } => {
                    for (offset, index) in input_of(*a).enumerate() {
                        drivers[index] =
                            Driver::Linear(vec![(output_of(*b).start + offset, *ratio)]);
                    }
                    for (offset, index) in input_of(*b).enumerate() {
                        drivers[index] =
                            Driver::Linear(vec![(output_of(*a).start + offset, -*ratio)]);
                    }
                }
                Connection::PowerMap { members, b } => {
                    let rows = self.spec(*b).width();
                    for MapMember { a, map, effort_map } in members {
                        let columns = self.spec(*a).width();
                        for (column, index) in input_of(*a).enumerate() {
                            let terms = (0..rows)
                                .filter_map(|row| {
                                    let coefficient = effort_map
                                        .as_ref()
                                        .map_or(map[row * columns + column], |adjoint| {
                                            adjoint[column * rows + row]
                                        });
                                    (coefficient != 0.0)
                                        .then_some((output_of(*b).start + row, coefficient))
                                })
                                .collect();
                            drivers[index] = Driver::Linear(terms);
                        }
                    }
                    // The members' terms follow one another in member order;
                    // one member gives the plain map's terms, in its order.
                    for (row, index) in input_of(*b).enumerate() {
                        let terms = members
                            .iter()
                            .flat_map(|MapMember { a, map, .. }| {
                                let columns = self.spec(*a).width();
                                let start = output_of(*a).start;
                                (0..columns).filter_map(move |column| {
                                    let coefficient = map[row * columns + column];
                                    (coefficient != 0.0).then_some((start + column, -coefficient))
                                })
                            })
                            .collect();
                        drivers[index] = Driver::Linear(terms);
                    }
                }
                Connection::HeatRoute {
                    source,
                    sinks,
                    remainder,
                    temperature,
                } => {
                    let source_out = output_of(*source).start;
                    for (sink, fraction) in sinks {
                        let share = (source_out, -*fraction);
                        match &mut drivers[input_of(*sink).start] {
                            // Another route already delivers here.
                            Driver::Linear(shares) => shares.push(share),
                            slot => *slot = Driver::Linear(vec![share]),
                        }
                    }
                    drivers[input_of(*source).start] = match temperature {
                        RouteTemperature::Fixed(kelvin) => Driver::Fixed(*kelvin),
                        RouteTemperature::SinkEquivalent => {
                            let terms: Vec<(usize, f64)> = sinks
                                .iter()
                                .map(|(sink, fraction)| (output_of(*sink).start, *fraction))
                                .collect();
                            let (fixed, extra) = remainder
                                .map_or((0.0, 0.0), |(_, fraction, kelvin)| {
                                    (fraction / kelvin, fraction)
                                });
                            let total = sinks.iter().map(|(_, f)| *f).sum::<f64>() + extra;
                            Driver::EquivalentTemperature {
                                terms,
                                fixed,
                                total,
                            }
                        }
                    };
                }
                Connection::Signal { from, to } => {
                    let from_out = output_of(*from).start;
                    for (offset, index) in input_of(*to).enumerate() {
                        drivers[index] = Driver::Linear(vec![(from_out + offset, 1.0)]);
                    }
                }
            }
        }
        for (boundary, (port, _)) in self.boundaries.iter().enumerate() {
            for (component, index) in input_of(*port).enumerate() {
                drivers[index] = Driver::Boundary {
                    boundary,
                    component,
                };
            }
        }

        // Which component owns each output and each input.
        let owner = |starts: &[usize], total: usize| -> Vec<usize> {
            let mut owners = Vec::with_capacity(total);
            for node in 0..count {
                owners.extend(core::iter::repeat_n(node, starts[node + 1] - starts[node]));
            }
            owners
        };
        let output_owner = owner(&output_start, outputs);
        let input_owner = owner(&input_start, inputs);

        // The end of a step at the ports of a load line. A port that names
        // an end input is driven there by the sum that drives its flow,
        // taken over what its neighbours return as the END of their flows:
        // the same outputs' companions, with the same signs and ratios. A
        // port that names none is joined only to neighbours whose flow ends
        // the step by the midpoint relation.
        for node in 0..count {
            let Some(line_ports) = self.nodes[node].line_ports() else {
                continue;
            };
            let ends = self.nodes[node].line_ends();
            for (position, &port) in line_ports.iter().enumerate() {
                let refuse = |reason: &'static str| AdmissionError::LineEnd {
                    port: self.label(PortRef { node, port }),
                    reason,
                };
                let flows = input_of(PortRef { node, port });
                let end_inputs = match ends.get(position).copied().flatten() {
                    Some(end_port) => {
                        let fits = self.nodes[node].ports().get(end_port).is_some_and(|spec| {
                            spec.carrier() == Carrier::Signal && spec.role() == PortRole::SignalIn
                        });
                        if !fits || self.used[node][end_port] {
                            return Err(refuse(
                                "its end input must be a signal input that nothing is connected \
                                 to: the kernel makes that wire",
                            ));
                        }
                        let range = input_of(PortRef {
                            node,
                            port: end_port,
                        });
                        if range.len() != flows.len() {
                            return Err(refuse("its end input has not the width of the port"));
                        }
                        Some(range)
                    }
                    None => None,
                };
                for (offset, flow) in flows.enumerate() {
                    let Driver::Linear(terms) = &drivers[flow] else {
                        // A flow prescribed from outside ends its step as
                        // it is held, or by the span the component is told.
                        if end_inputs.is_some() {
                            return Err(refuse(
                                "its flow is prescribed, and a prescribed flow has no end to \
                                 tell: the port takes no end input",
                            ));
                        }
                        continue;
                    };
                    let terms = terms.clone();
                    let mut end_terms = Vec::with_capacity(terms.len());
                    let (mut told, mut midpoint) = (0_usize, 0_usize);
                    for (source, coefficient) in terms {
                        let from = output_owner[source];
                        let local = source - output_start[from];
                        let from_port = (0..self.nodes[from].ports().len())
                            .find(|&port| self.layouts[from].output_range(port).contains(&local))
                            .expect("every output belongs to a port of its component");
                        let within = local - self.layouts[from].output_range(from_port).start;
                        match self.nodes[from].end_of_step(from_port) {
                            EndOfStep::Midpoint => midpoint += 1,
                            EndOfStep::Untold => {
                                return Err(refuse(
                                    "a neighbour's flow does not end its step by the midpoint \
                                     relation, and it does not tell its end",
                                ));
                            }
                            EndOfStep::Told(companion) => {
                                let fits =
                                    self.nodes[from].ports().get(companion).is_some_and(|spec| {
                                        spec.carrier() == Carrier::Signal
                                            && spec.role() == PortRole::SignalOut
                                    });
                                let range = if fits {
                                    self.layouts[from].output_range(companion)
                                } else {
                                    0..0
                                };
                                if range.len() != self.layouts[from].output_range(from_port).len() {
                                    return Err(refuse(
                                        "a neighbour names, for the end of its flow, a port \
                                         that is not a signal output of the flow's width",
                                    ));
                                }
                                told += 1;
                                end_terms
                                    .push((output_start[from] + range.start + within, coefficient));
                            }
                        }
                    }
                    match (&end_inputs, told, midpoint) {
                        (None, 0, _) => {}
                        (None, _, _) => {
                            return Err(refuse(
                                "a neighbour tells the end of its flow, which does not end its \
                                 step by the midpoint relation, and the port takes no end input",
                            ));
                        }
                        (Some(_), 0, _) => {
                            return Err(refuse(
                                "it takes an end input, and no neighbour tells the end of its \
                                 flow: they end their steps by the midpoint relation",
                            ));
                        }
                        (Some(range), _, 0) => {
                            drivers[range.start + offset] = Driver::Linear(end_terms);
                        }
                        (Some(_), _, _) => {
                            return Err(refuse(
                                "some of its neighbours tell the end of their flow and others \
                                 do not",
                            ));
                        }
                    }
                }
            }
        }

        let mut successors: Vec<Vec<usize>> = vec![Vec::new(); count];
        let mut sources = Vec::new();
        for (index, driver) in drivers.iter().enumerate() {
            sources.clear();
            driver.sources(&mut sources);
            for &source in &sources {
                let (from, to) = (output_owner[source], input_owner[index]);
                if !successors[from].contains(&to) {
                    successors[from].push(to);
                }
            }
        }
        for list in &mut successors {
            list.sort_unstable();
        }
        let mut block_of = vec![0_usize; count];
        let node_blocks = blocks_in_dependency_order(&successors);
        for (block, nodes) in node_blocks.iter().enumerate() {
            for &node in nodes {
                block_of[node] = block;
            }
        }
        let mut blocks: Vec<Block> = node_blocks
            .into_iter()
            .map(|nodes| Block {
                nodes,
                unknowns: Vec::new(),
                known: Vec::new(),
            })
            .collect();
        for (index, driver) in drivers.iter().enumerate() {
            sources.clear();
            driver.sources(&mut sources);
            let block = block_of[input_owner[index]];
            if sources
                .iter()
                .any(|&source| block_of[output_owner[source]] == block)
            {
                blocks[block].unknowns.push(index);
            }
        }
        for block in &mut blocks {
            for &node in &block.nodes {
                for index in input_start[node]..input_start[node + 1] {
                    if !block.unknowns.contains(&index) {
                        block.known.push(index);
                    }
                }
            }
        }
        // A component on a load line whose flows are solved in a loop: the
        // drivers of those flows gain a term per effort of the component,
        // for the line's compliance, and the block keeps the line.
        let mut lines: Vec<Option<Line>> = Vec::new();
        lines.resize_with(blocks.len(), || None);
        let mut line_of = vec![usize::MAX; count];
        for node in 0..count {
            let Some(ports) = self.nodes[node].line_ports() else {
                continue;
            };
            let block = block_of[node];
            let mut places = Vec::new();
            for &port in ports {
                let flows = self.layouts[node].input_range(port);
                let efforts = self.layouts[node].output_range(port);
                for (flow, effort) in flows.zip(efforts) {
                    let input = input_start[node] + flow;
                    let unknown = blocks[block]
                        .unknowns
                        .iter()
                        .position(|&index| index == input)
                        .unwrap_or(usize::MAX);
                    places.push((input, output_start[node] + effort, unknown));
                }
            }
            if places.iter().all(|&(_, _, unknown)| unknown == usize::MAX) {
                continue;
            }
            if let Some(first) = &lines[block] {
                return Err(AdmissionError::TwoLines {
                    first: self.names[first.node].clone(),
                    second: self.names[node].clone(),
                });
            }
            let own = places.len();
            for &(input, _, unknown) in &places {
                if unknown == usize::MAX {
                    continue;
                }
                if let Driver::Linear(terms) = &mut drivers[input] {
                    let mut terms = core::mem::take(terms);
                    terms.extend(places.iter().map(|&(_, effort, _)| (effort, 0.0)));
                    drivers[input] = Driver::OnLine { terms, own };
                }
            }
            // The end inputs wired above are unknowns of the block like
            // the flows, and the line has a row for each.
            let mut ends = Vec::new();
            for (position, _) in ports.iter().enumerate() {
                let Some(end_port) = self.nodes[node]
                    .line_ends()
                    .get(position)
                    .copied()
                    .flatten()
                else {
                    continue;
                };
                for end in self.layouts[node].input_range(end_port) {
                    let input = input_start[node] + end;
                    let unknown = blocks[block]
                        .unknowns
                        .iter()
                        .position(|&index| index == input)
                        .unwrap_or(usize::MAX);
                    if unknown != usize::MAX
                        && let Driver::Linear(terms) = &mut drivers[input]
                    {
                        let mut terms = core::mem::take(terms);
                        terms.extend(places.iter().map(|&(_, effort, _)| (effort, 0.0)));
                        drivers[input] = Driver::OnLine { terms, own };
                    }
                    ends.push((input, unknown));
                }
            }
            line_of[node] = block;
            lines[block] = Some(Line::new(node, places, ends));
        }
        let block_work = blocks
            .iter()
            .zip(&lines)
            .map(|(block, line)| {
                // A part without a declared tangent is differenced, and a
                // difference quotient is rounded at the point it is taken
                // at: its slopes are not the same numbers from one step to
                // the next, so factors made from them are not kept. (A
                // component on a load line has no declared tangent, so a
                // block with one keeps none either.)
                let affine = block.nodes.iter().all(|&node| {
                    self.nodes[node].class() == Class::ExactLinear && self.nodes[node].has_tangent()
                });
                let mut work =
                    BlockWork::plan(block, &drivers, &input_owner, &output_start, affine);
                if line.is_some() {
                    let declares = |node: usize| self.nodes[node].has_tangent();
                    work.found =
                        Found::plan(block, &input_owner, &input_start, &output_start, &declares);
                }
                work
            })
            .collect();

        // Tolerance group and label of every input. Power and heat inputs of
        // one carrier and side are compared with the largest of their kind in
        // the block. The coefficients of a stream port are different
        // quantities (a pressure, an enthalpy, a mass fraction), so each
        // position is a kind of its own. A signal has no kind, so each scalar
        // signal input is a group of its own and is judged against its own
        // size only.
        let mut groups = Vec::new();
        let mut signal_inputs = Vec::new();
        let mut input_group = vec![0_usize; inputs];
        let mut input_label = vec![String::new(); inputs];
        let mut input_is_temperature = vec![false; inputs];
        let mut input_is_energy_rate = vec![false; inputs];
        for node in 0..count {
            for (port, spec) in self.nodes[node].ports().iter().enumerate() {
                let reference = PortRef { node, port };
                for (offset, index) in input_of(reference).enumerate() {
                    if spec.carrier() == Carrier::Signal {
                        signal_inputs.push(index);
                    } else {
                        input_group[index] = tolerance_group(&mut groups, spec, offset);
                    }
                    input_label[index] = self.label(reference);
                    input_is_temperature[index] =
                        spec.carrier() == Carrier::Heat && spec.role() == PortRole::EffortIn;
                    input_is_energy_rate[index] = match (spec.carrier(), spec.role()) {
                        (Carrier::Heat, PortRole::FlowIn) => true,
                        (Carrier::Stream, PortRole::FlowIn) => offset == STREAM_ENTHALPY_RATE,
                        _ => false,
                    };
                }
            }
        }
        let mut group_count = groups.len();
        for index in signal_inputs {
            input_group[index] = group_count;
            group_count += 1;
        }
        // Every claimed port got a driver above; an input without one would
        // be read as NaN by the first step, so it is refused here instead.
        if let Some(index) = drivers
            .iter()
            .position(|driver| matches!(driver, Driver::Unset))
        {
            return Err(AdmissionError::Unconnected {
                port: input_label[index].clone(),
            });
        }

        let ledger = Ledger {
            nodes: self
                .names
                .iter()
                .enumerate()
                .map(|(index, name)| NodeAccount {
                    name: name.clone(),
                    carried: vec![None; self.nodes[index].ports().iter()
                        .map(PortSpec::carried_quantities).max().unwrap_or(0)],
                    has_streams: self.nodes[index]
                        .ports()
                        .iter()
                        .any(|port| port.carrier() == Carrier::Stream),
                    ..NodeAccount::default()
                })
                .collect(),
            connections: self.connections.iter().map(|connection| {
                let mut members = Vec::new();
                connection.energy_members(&mut members);
                let carried = members.iter().map(|port| self.spec(*port).carried_quantities())
                    .max().unwrap_or(0);
                ConnectionAccount {
                    carried_defect: vec![Sum::default(); carried],
                    ..ConnectionAccount::default()
                }
            }).collect(),
            boundaries: self
                .boundaries
                .iter()
                .map(|(port, _)| BoundaryAccount {
                    port: self.label(*port),
                    carried_in: vec![Sum::default(); self.spec(*port).carried_quantities()],
                    ..BoundaryAccount::default()
                })
                .collect(),
            sinks: self
                .sinks
                .iter()
                .map(|sink| SinkAccount {
                    name: sink.name.as_str().to_owned(),
                    ..SinkAccount::default()
                })
                .collect(),
            ..Ledger::default()
        };
        let boundary_values = self
            .boundaries
            .iter()
            .map(|(port, _)| vec![0.0; self.spec(*port).input_width()])
            .collect();
        let widest_output = self
            .layouts
            .iter()
            .map(PortLayout::output_len)
            .max()
            .unwrap_or(0);
        let widest_input = self
            .layouts
            .iter()
            .map(PortLayout::input_len)
            .max()
            .unwrap_or(0);
        // Heat routes that share a receiver are audited together: that
        // port's energy is the sum of what each of them delivered. The group
        // is kept under its first connection; the others keep an empty list.
        let mut audit_groups: Vec<Vec<usize>> = (0..self.connections.len())
            .map(|index| vec![index])
            .collect();
        let mut leader: Vec<usize> = (0..self.connections.len()).collect();
        for later in 0..self.connections.len() {
            let Connection::HeatRoute { sinks: mine, .. } = &self.connections[later] else {
                continue;
            };
            for earlier in 0..later {
                let Connection::HeatRoute { sinks: theirs, .. } = &self.connections[earlier] else {
                    continue;
                };
                let shares = mine
                    .iter()
                    .any(|(port, _)| theirs.iter().any(|(other, _)| other == port));
                let (keep, merge) = (
                    leader[earlier].min(leader[later]),
                    leader[earlier].max(leader[later]),
                );
                if shares && keep != merge {
                    let moved = core::mem::take(&mut audit_groups[merge]);
                    for &index in &moved {
                        leader[index] = keep;
                    }
                    audit_groups[keep].extend(moved);
                    audit_groups[keep].sort_unstable();
                }
            }
        }

        let port_carried = self
            .nodes
            .iter()
            .flat_map(|node| {
                node.ports()
                    .iter()
                    .map(|port| vec![0.0; port.carried_quantities()])
            })
            .collect();
        let carried_audits = ledger
            .nodes
            .iter()
            .map(|node| vec![None; node.carried.len()])
            .collect();
        let connection_carried = ledger
            .connections
            .iter()
            .map(|connection| vec![0.0; connection.carried_defect.len()])
            .collect();
        Ok(System {
            identity: Rc::new(()),
            names: self.names,
            nodes: self.nodes,
            layouts: self.layouts,
            accepted: self.accepted,
            candidate: self.candidate,
            connections: self.connections,
            audit_groups,
            boundaries: self.boundaries,
            boundary_values,
            drivers,
            blocks,
            block_work,
            lines,
            line_of,
            step_scratch: StepScratch::default(),
            input_start,
            output_start,
            port_start,
            input_owner,
            input_group,
            input_label,
            input_is_temperature,
            input_is_energy_rate,
            group_count,
            inputs: vec![0.0; inputs],
            outputs: vec![0.0; outputs],
            accepted_inputs: vec![0.0; inputs],
            accepted_outputs: vec![0.0; outputs],
            input_allowance: vec![0.0; inputs],
            accepted_allowance: vec![0.0; inputs],
            receipts: vec![Receipt::lossless(0.0, 0.0); count],
            port_energy: vec![0.0; ports],
            port_entropy: vec![0.0; ports],
            port_mass: vec![0.0; ports],
            port_carried,
            carried_audits,
            connection_carried,
            scratch_outputs: vec![0.0; widest_output],
            scratch_direction: vec![0.0; widest_input],
            scratch_first: vec![0.0; widest_output],
            scratch_kept: vec![0.0; widest_output],
            scratch_entries: vec![Entry::First; widest_output],
            unresolved: 0,
            time: Sum::default(),
            steps: 0,
            ledger,
            config: SolveConfig::default(),
        })
    }
}

enum Trial {
    Done { iterations: u32 },
    Shorten { dt_s: f64 },
}

/// Lists a step fills and empties, kept between steps.
#[derive(Default)]
struct StepScratch {
    defects: Vec<f64>,
    mass_audits: Vec<Option<MassStepAudit>>,
    members: Vec<PortRef>,
    group_members: Vec<PortRef>,
    connection_rows: Vec<(f64, f64, f64)>,
    sink_heat: Vec<f64>,
}

/// What a looped block keeps between steps: the workspaces of its solve,
/// which residual rows read which component's outputs, and, for a block
/// that is affine throughout, the factors of its interface Jacobian.
#[derive(Default)]
struct BlockWork {
    residual: Vec<f64>,
    /// First-order size of the rounding in each residual: how strongly its
    /// driver's sources depend on the block's inputs, times the size of
    /// those inputs. Known after a tangent assembly.
    rounding: Vec<f64>,
    allowance: Vec<f64>,
    group_scale: Vec<f64>,
    gradient: Vec<(usize, f64)>,
    /// The interface Jacobian and its factors: the unknowns that are
    /// copies substituted, the rest eliminated together.
    condensed: Condensed,
    /// Per unknown, the part of `readers` that reads its owner's outputs.
    column_readers: Vec<Range<usize>>,
    /// `(row, position in the owner's outputs, slope)` for every term of a
    /// driver of the block, grouped by the component whose output it
    /// reads, in row order and then in the driver's order. The slope of a
    /// driver that is not linear is remade before every assembly.
    readers: Vec<(usize, usize, f64)>,
    /// Per unknown, where the Jacobian keeps the entry of each of its
    /// readers: `slots[slot_start[column]..]`, in the readers' order.
    slots: Vec<usize>,
    slot_start: Vec<usize>,
    /// Rows whose driver's slopes depend on the outputs it reads, each with
    /// its part of `varying_terms`: per term of the driver, the reader that
    /// carries its slope (`usize::MAX` for a term that reads no output of
    /// this block; a heat route's receivers are always in its source's).
    varying: Vec<(usize, Range<usize>)>,
    varying_terms: Vec<usize>,
    /// Every component is exactly linear and declares its tangent, and
    /// every driver is linear, so the Jacobian is taken to depend on the
    /// step length and on nothing else: the kept factors and slopes are
    /// then the numbers a fresh assembly would form, and using them does
    /// not let an earlier step decide a later one.
    affine: bool,
    /// Of an affine block: the step its factors were made for (NaN if it
    /// has none).
    kept_dt_s: f64,
    /// The Jacobian's entries off the identity as `(row, column, value)`,
    /// in the order the most recent assembly formed them: how the target
    /// of wire `row` moves with input `column`. An affine block reuses
    /// them with its factors at an unchanged step length. The second
    /// convergence test reads them, and only once the step it is judging
    /// has assembled or adopted them: nothing else crosses a step.
    slopes: Vec<(usize, usize, f64)>,
    /// Per wire, for the second convergence test: the rounding floor of
    /// its residual at the current inputs, and what the residuals of the
    /// other wires make of its target when carried once through `slopes`.
    floor: Vec<f64>,
    carried: Vec<f64>,
    /// Per wire, the size its residual is relative to: that of its
    /// tolerance group, or its own in a block with a load line.
    scale: Vec<f64>,
    /// Per wire of a block with a load line that is not a flow on the
    /// line: the rounding of the largest wire of its kind in the block.
    /// Zero for every other wire.
    kind_rounding: Vec<f64>,
    /// What finding the block's load line evaluated, kept for the solve's
    /// first pass. Empty for a block without a load line.
    found: Found,
}

impl BlockWork {
    /// Workspaces and reader lists of `block`. `affine` says whether every
    /// component of it is exactly linear and declares its tangent.
    fn plan(
        block: &Block,
        drivers: &[Driver],
        input_owner: &[usize],
        output_start: &[usize],
        affine: bool,
    ) -> Self {
        let m = block.unknowns.len();
        // Rows whose driver is not linear, and room for where its terms go.
        let mut varying = Vec::new();
        let mut varying_terms = Vec::new();
        let mut terms_of_row = vec![usize::MAX; m];
        let mut sources = Vec::new();
        for (row, &target) in block.unknowns.iter().enumerate() {
            if matches!(drivers[target], Driver::Linear(_)) {
                continue;
            }
            sources.clear();
            drivers[target].sources(&mut sources);
            terms_of_row[row] = varying_terms.len();
            varying_terms.resize(varying_terms.len() + sources.len(), usize::MAX);
            varying.push((row, terms_of_row[row]..varying_terms.len()));
        }
        let mut readers = Vec::new();
        let mut of_node = Vec::with_capacity(block.nodes.len());
        for &node in &block.nodes {
            let outputs = output_start[node]..output_start[node + 1];
            let start = readers.len();
            for (row, &target) in block.unknowns.iter().enumerate() {
                if let Driver::Linear(terms) = &drivers[target] {
                    for &(source, slope) in terms {
                        if outputs.contains(&source) {
                            readers.push((row, source - outputs.start, slope));
                        }
                    }
                } else {
                    sources.clear();
                    drivers[target].sources(&mut sources);
                    for (term, &source) in sources.iter().enumerate() {
                        if outputs.contains(&source) {
                            varying_terms[terms_of_row[row] + term] = readers.len();
                            readers.push((row, source - outputs.start, 0.0));
                        }
                    }
                }
            }
            of_node.push(start..readers.len());
        }
        let column_readers: Vec<Range<usize>> = block
            .unknowns
            .iter()
            .map(|&index| {
                let position = block
                    .nodes
                    .binary_search(&input_owner[index])
                    .expect("an unknown belongs to a component of its block");
                of_node[position].clone()
            })
            .collect();
        // Column `j` of the Jacobian has an entry in every row that reads
        // an output of the component unknown `j` is an input of.
        let pattern: Vec<Vec<usize>> = column_readers
            .iter()
            .map(|range| {
                let mut rows: Vec<usize> = readers[range.clone()]
                    .iter()
                    .map(|&(row, _, _)| row)
                    .collect();
                rows.sort_unstable();
                rows.dedup();
                rows
            })
            .collect();
        let condensed = Condensed::plan(&pattern);
        let mut slots = Vec::new();
        let mut slot_start = Vec::with_capacity(m);
        for (column, range) in column_readers.iter().enumerate() {
            slot_start.push(slots.len());
            slots.extend(
                readers[range.clone()]
                    .iter()
                    .map(|&(row, _, _)| condensed.slot(row, column)),
            );
        }
        let entries = slots.len();
        Self {
            residual: vec![0.0; m],
            rounding: vec![0.0; m],
            allowance: vec![0.0; m],
            group_scale: Vec::new(),
            gradient: Vec::new(),
            condensed,
            column_readers,
            readers,
            slots,
            slot_start,
            affine: affine && varying.is_empty(),
            varying,
            varying_terms,
            kept_dt_s: f64::NAN,
            // One entry per Jacobian entry: an assembly never grows it.
            slopes: Vec::with_capacity(entries),
            floor: vec![0.0; m],
            carried: vec![0.0; m],
            scale: vec![0.0; m],
            kind_rounding: vec![0.0; m],
            found: Found::default(),
        }
    }
}

/// What finding a block's load line evaluated and still holds at the solve's
/// first pass (`solve_lent`): per component of the block that declares its
/// tangent, the inputs finding the line proposed it at, and per unknown the
/// column that tangent gave. Sized at seal, never in a step.
#[derive(Default)]
struct Found {
    /// Per position in the block's components: whether it declares its
    /// tangent (only such a component is kept), and whether it was kept.
    declares: Vec<bool>,
    proposed: Vec<bool>,
    /// Where each position's inputs start in `inputs`.
    input_start: Vec<usize>,
    inputs: Vec<f64>,
    /// Per unknown: the position of its owner in the block's components,
    /// whether its column was kept, and where it starts in `columns`.
    owner_position: Vec<usize>,
    column: Vec<bool>,
    column_start: Vec<usize>,
    columns: Vec<f64>,
    /// Per position, in the first pass of a step: kept at bit for bit the
    /// inputs it holds now, so its proposal and its columns are those.
    repeated: Vec<bool>,
}

impl Found {
    /// Room for every component's inputs and every unknown's column (its
    /// owner's outputs) of `block`; `declares` says whether a component
    /// declares its tangent.
    fn plan(
        block: &Block,
        input_owner: &[usize],
        input_start: &[usize],
        output_start: &[usize],
        declares: &dyn Fn(usize) -> bool,
    ) -> Self {
        let (mut starts, mut inputs) = (Vec::new(), 0_usize);
        for &node in &block.nodes {
            starts.push(inputs);
            inputs += input_start[node + 1] - input_start[node];
        }
        let (mut owner_position, mut column_start, mut columns) = (Vec::new(), Vec::new(), 0_usize);
        for &index in &block.unknowns {
            let owner = input_owner[index];
            owner_position.push(
                block
                    .nodes
                    .binary_search(&owner)
                    .expect("an unknown belongs to a component of its block"),
            );
            column_start.push(columns);
            columns += output_start[owner + 1] - output_start[owner];
        }
        let (components, unknowns) = (block.nodes.len(), block.unknowns.len());
        Self {
            declares: block.nodes.iter().map(|&node| declares(node)).collect(),
            proposed: vec![false; components],
            input_start: starts,
            inputs: vec![0.0; inputs],
            owner_position,
            column: vec![false; unknowns],
            column_start,
            columns: vec![0.0; columns],
            repeated: vec![false; components],
        }
    }

    /// Forget what an earlier finding of the line kept.
    fn clear(&mut self) {
        self.proposed.fill(false);
        self.column.fill(false);
    }

    /// Keep `inputs`, at which the component at `position` was just
    /// proposed, if it declares its tangent.
    fn keep_inputs(&mut self, position: usize, inputs: &[f64]) {
        if self.declares[position] {
            let start = self.input_start[position];
            self.inputs[start..start + inputs.len()].copy_from_slice(inputs);
            self.proposed[position] = true;
        }
    }

    /// Keep `values`, the column of unknown `column` its owner's tangent
    /// gave, if the owner declares its tangent (a differenced column is not
    /// kept: differencing writes the owner's candidate).
    fn keep_column(&mut self, column: usize, values: &[f64]) {
        if self.declares[self.owner_position[column]] {
            let start = self.column_start[column];
            self.columns[start..start + values.len()].copy_from_slice(values);
            self.column[column] = true;
        }
    }

    /// Whether the component at `position` was kept at bit for bit `inputs`;
    /// recorded in `repeated` for the first assembly.
    fn repeats(&mut self, position: usize, inputs: &[f64]) -> bool {
        let start = self.input_start[position];
        self.repeated[position] = self.proposed[position]
            && inputs
                .iter()
                .zip(&self.inputs[start..start + inputs.len()])
                .all(|(now, then)| now.to_bits() == then.to_bits());
        self.repeated[position]
    }

    /// The kept column of unknown `column`, `width` long, if its owner
    /// repeats in this first pass.
    fn column_of(&self, column: usize, width: usize) -> Option<&[f64]> {
        if self.repeated[self.owner_position[column]] && self.column[column] {
            let start = self.column_start[column];
            Some(&self.columns[start..start + width])
        } else {
            None
        }
    }
}

/// Components joined by connections, stepped as one.
///
/// The system owns every state in two generations, accepted and candidate. A
/// step fills the candidates, audits and admits them, and only then publishes
/// all of them by swapping; a failed step leaves every accepted state, the
/// clock and the ledger untouched.
pub struct System {
    /// Allocation identity survives moves and cannot be reused while a
    /// checkpoint exists. Equal graph descriptions still get distinct tokens.
    identity: Rc<()>,
    names: Vec<String>,
    nodes: Vec<Box<dyn Node>>,
    layouts: Vec<PortLayout>,
    accepted: Vec<Box<dyn Any>>,
    candidate: Vec<Box<dyn Any>>,
    connections: Vec<Connection>,
    /// Per connection, the connections audited in its row (itself, plus the
    /// heat routes that share a receiver with it); empty if an earlier
    /// connection's row holds it.
    audit_groups: Vec<Vec<usize>>,
    boundaries: Vec<(PortRef, BoundaryFn)>,
    boundary_values: Vec<Vec<f64>>,
    drivers: Vec<Driver>,
    blocks: Vec<Block>,
    /// Per block, what its solve keeps between steps.
    block_work: Vec<BlockWork>,
    /// Per block, its component on a load line, if its flows are solved
    /// there; and per component, the block whose line it is on
    /// (`usize::MAX` for none).
    lines: Vec<Option<Line>>,
    line_of: Vec<usize>,
    step_scratch: StepScratch,
    input_start: Vec<usize>,
    output_start: Vec<usize>,
    port_start: Vec<usize>,
    input_owner: Vec<usize>,
    input_group: Vec<usize>,
    input_label: Vec<String>,
    input_is_temperature: Vec<bool>,
    /// Inputs that are a heat rate or a stream's enthalpy rate into a
    /// storage: energy per second, made by other components from products
    /// of their efforts and flows.
    input_is_energy_rate: Vec<bool>,
    group_count: usize,
    inputs: Vec<f64>,
    outputs: Vec<f64>,
    accepted_inputs: Vec<f64>,
    accepted_outputs: Vec<f64>,
    /// Per scalar input, the part of the residual the block solve accepted
    /// for it in the current trial that no relative tolerance covers: its
    /// rounding estimate and the configured absolute tolerance, and in a
    /// block with a load line the rounding of the largest wire of its kind
    /// in the block or, for a flow on the line, of the difference it is.
    /// Zero for an input copied from its driver.
    input_allowance: Vec<f64>,
    /// `input_allowance` of the last accepted step; zero after a restore.
    accepted_allowance: Vec<f64>,
    receipts: Vec<Receipt>,
    port_energy: Vec<f64>,
    port_entropy: Vec<f64>,
    port_mass: Vec<f64>,
    /// Carried-account workspaces are allocated at seal, never per step.
    port_carried: Vec<Vec<f64>>,
    carried_audits: Vec<Vec<Option<CarriedStepAudit>>>,
    connection_carried: Vec<Vec<f64>>,
    scratch_outputs: Vec<f64>,
    scratch_direction: Vec<f64>,
    /// A differenced column as the first step left it, as its entries are
    /// resolved, and where each stands.
    scratch_first: Vec<f64>,
    scratch_kept: Vec<f64>,
    scratch_entries: Vec<Entry>,
    /// Entries of differenced columns left unresolved in this step.
    unresolved: u32,
    time: Sum,
    steps: u64,
    ledger: Ledger,
    config: SolveConfig,
}

impl System {
    /// Capture the accepted state for later inspection or rewind.
    ///
    /// Capture also works after a refused or cancelled step: trial candidates
    /// and scratch values are excluded. The system itself is unchanged.
    ///
    /// # Errors
    /// A component state has the wrong erased type (a kernel invariant breach).
    pub fn checkpoint(&self) -> Result<Checkpoint, CheckpointError> {
        Ok(Checkpoint {
            identity: Rc::clone(&self.identity),
            accepted: self.clone_states(&self.accepted)?,
            inputs: self.accepted_inputs.clone(),
            outputs: self.accepted_outputs.clone(),
            time: self.time,
            steps: self.steps,
            ledger: self.ledger.clone(),
            config: self.config,
        })
    }

    /// Restore an accepted-state checkpoint without consuming it.
    ///
    /// Every saved state is cloned before anything is published, so a refusal
    /// leaves the system untouched. Candidate slots remain disposable scratch;
    /// the next proposal overwrites them under the Component contract.
    /// The graph, boundary closures and external execution context are retained.
    ///
    /// # Errors
    /// A checkpoint from another sealed system, or an erased state type mismatch.
    pub fn restore(&mut self, checkpoint: &Checkpoint) -> Result<(), CheckpointError> {
        if !Rc::ptr_eq(&self.identity, &checkpoint.identity) {
            return Err(CheckpointError::ForeignSystem);
        }
        let accepted = self.clone_states(&checkpoint.accepted)?;
        let inputs = checkpoint.inputs.clone();
        let outputs = checkpoint.outputs.clone();
        let ledger = checkpoint.ledger.clone();

        self.accepted = accepted;
        self.accepted_inputs = inputs;
        self.accepted_outputs = outputs;
        self.accepted_allowance.fill(0.0);
        self.inputs.copy_from_slice(&self.accepted_inputs);
        self.outputs.copy_from_slice(&self.accepted_outputs);
        self.time = checkpoint.time;
        self.steps = checkpoint.steps;
        self.ledger = ledger;
        self.config = checkpoint.config;
        Ok(())
    }

    fn clone_states(&self, states: &[Box<dyn Any>]) -> Result<Vec<Box<dyn Any>>, CheckpointError> {
        self.nodes
            .iter()
            .zip(states)
            .enumerate()
            .map(|(index, (node, state))| {
                node.clone_state(state.as_ref())
                    .map_err(|source| CheckpointError::State {
                        node: self.names[index].clone(),
                        source,
                    })
            })
            .collect()
    }

    /// Replace the driver's limits.
    pub fn set_config(&mut self, config: SolveConfig) {
        self.config = config;
    }

    /// Simulated time [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.time.value()
    }

    /// Accepted steps so far.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// The window ledger.
    #[must_use]
    pub fn ledger(&self) -> &Ledger {
        &self.ledger
    }

    /// Number of blocks in the plan, and how many of them are loops.
    #[must_use]
    pub fn plan_shape(&self) -> (usize, usize) {
        let loops = self.blocks.iter().filter(|block| block.looped()).count();
        (self.blocks.len(), loops)
    }

    /// Per loop of the plan, in its order: how many inputs it solves for,
    /// and how many of those are solved together. The others are copies and
    /// follow from them by substitution.
    #[must_use]
    pub fn plan_loops(&self) -> Vec<(usize, usize)> {
        self.blocks
            .iter()
            .zip(&self.block_work)
            .filter(|(block, _)| block.looped())
            .map(|(block, work)| (block.unknowns.len(), work.condensed.kept()))
            .collect()
    }

    /// The accepted state of `node`, if it has type `S`.
    #[must_use]
    pub fn state<S: 'static>(&self, node: NodeId) -> Option<&S> {
        self.accepted.get(node.0)?.downcast_ref::<S>()
    }

    /// The port of `node` with identifier `id`.
    #[must_use]
    pub fn port(&self, node: NodeId, id: &str) -> Option<PortRef> {
        let position = self
            .nodes
            .get(node.0)?
            .ports()
            .iter()
            .position(|port| port.id().as_str() == id)?;
        Some(PortRef {
            node: node.0,
            port: position,
        })
    }

    fn input_range(&self, port: PortRef) -> Range<usize> {
        let range = self.layouts[port.node].input_range(port.port);
        self.input_start[port.node] + range.start..self.input_start[port.node] + range.end
    }

    fn output_range(&self, port: PortRef) -> Range<usize> {
        let range = self.layouts[port.node].output_range(port.port);
        self.output_start[port.node] + range.start..self.output_start[port.node] + range.end
    }

    /// What `port` received during the last accepted step.
    #[must_use]
    pub fn last_input(&self, port: PortRef) -> &[f64] {
        &self.accepted_inputs[self.input_range(port)]
    }

    /// What `port` returned during the last accepted step.
    #[must_use]
    pub fn last_output(&self, port: PortRef) -> &[f64] {
        &self.accepted_outputs[self.output_range(port)]
    }

    /// What the solve of the last accepted step told the interface audit
    /// it had left in each input of `port`, in the input's unit: its
    /// rounding estimate and the configured absolute tolerance; in a block
    /// with a load line also the rounding of the largest wire of its kind
    /// in the block, or for a flow on the line the rounding of the
    /// difference it is. Never the solve's relative tolerance. Zero for an
    /// input that no block solved for, and after [`Self::restore`] until a
    /// step is taken: a checkpoint does not carry it.
    #[must_use]
    pub fn last_allowance(&self, port: PortRef) -> &[f64] {
        &self.accepted_allowance[self.input_range(port)]
    }

    /// Starting value for the input side of `port` in the next step's
    /// interface solve. Temperatures inside a loop need one before the first
    /// step. Returns `false`, and changes nothing, if `values` does not have
    /// the port's input width.
    pub fn guess(&mut self, port: PortRef, values: &[f64]) -> bool {
        let range = self.input_range(port);
        let fits = values.len() == range.len();
        if fits {
            self.accepted_inputs[range].copy_from_slice(values);
        }
        fits
    }

    fn node_inputs(&self, node: usize) -> Range<usize> {
        self.input_start[node]..self.input_start[node + 1]
    }

    fn node_outputs(&self, node: usize) -> Range<usize> {
        self.output_start[node]..self.output_start[node + 1]
    }

    fn names_of(&self, nodes: &[usize]) -> Vec<String> {
        nodes.iter().map(|&node| self.names[node].clone()).collect()
    }

    /// A component's refusal during a trial. Once cancellation has been
    /// requested every refusal is reported as [`SystemError::Cancelled`], so
    /// callers see one cancellation outcome whichever component noticed first.
    fn refused(&self, node: usize, error: OwnerError, cx: &Cx<'_>) -> SystemError {
        if cx.is_cancel_requested() {
            SystemError::Cancelled
        } else {
            SystemError::Component {
                node: self.names[node].clone(),
                phase: Phase::Propose,
                error,
            }
        }
    }

    /// Admit a component's request for a shorter step, or name the component
    /// that asked for a step that is not shorter.
    fn shorter(&self, node: usize, dt_s: f64, shorter: f64) -> Result<Option<f64>, SystemError> {
        if shorter > 0.0 && shorter < dt_s {
            Ok(Some(shorter))
        } else {
            Err(SystemError::Contract {
                node: self.names[node].clone(),
                violation: ContractViolation::Step { dt_s: shorter },
            })
        }
    }

    /// Propose `node` at the current inputs. `Ok(Some(dt))` asks to shorten.
    fn propose(&mut self, node: usize, dt_s: f64, cx: &Cx<'_>) -> Result<Option<f64>, SystemError> {
        let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
        // A component on a load line is stepped there: its flow inputs
        // hold the released flows of this step's line.
        let line = self.lines.get(self.line_of[node]).and_then(Option::as_ref);
        let proposal = if let Some(line) = line {
            self.nodes[node].propose_on_line(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs],
                &line.compliance,
                self.candidate[node].as_mut(),
                &mut self.outputs[outputs],
                cx,
            )
        } else {
            self.nodes[node].propose(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs],
                self.candidate[node].as_mut(),
                &mut self.outputs[outputs],
                cx,
            )
        };
        match proposal {
            Ok(Proposal::Step(receipt)) => {
                self.receipts[node] = receipt;
                Ok(None)
            }
            Ok(Proposal::Shorten { dt_s: shorter }) => self.shorter(node, dt_s, shorter),
            Err(error) => Err(self.refused(node, error, cx)),
        }
    }

    /// Column `local` of the port-map tangent of `node`, into
    /// `scratch_outputs`. `Ok(Some(dt))` asks to shorten. A differenced
    /// column is resolved against the rounding of its outputs
    /// (`resolve_column`), a load line's and a block solve's alike: the
    /// line component's relation is solved on the one, and an update is
    /// made with the other. An update made with an entry the rounding took
    /// moves the block by a fixed-point iteration in that entry's place,
    /// which need not contract: two capacities near rest behind a
    /// conductance, differenced, were refused at the iteration limit
    /// (`tests/elements.rs`, the capacity near rest), as were three cells of
    /// a timeline near rest with entries of 0 to 2 units.
    fn tangent_column(
        &mut self,
        node: usize,
        local: usize,
        dt_s: f64,
        group_scale: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SystemError> {
        let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
        let width = outputs.len();
        if self.nodes[node].has_tangent() {
            let direction = &mut self.scratch_direction[..inputs.len()];
            direction.fill(0.0);
            direction[local] = 1.0;
            let tangent = self.nodes[node].jvp(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs.clone()],
                &self.scratch_direction[..inputs.len()],
                &mut self.scratch_outputs[..width],
                cx,
            );
            return match tangent {
                Ok(()) => Ok(None),
                Err(error) => Err(self.refused(node, error, cx)),
            };
        }
        // A component added on a load line is a function of what it is
        // given, whatever its class, and is differenced as one: of the
        // released flows where the block solves for its line, and of its
        // other inputs where the line's flows are all given from outside
        // the block (a slider between prescribed velocities whose heat
        // warms a body).
        if self.nodes[node].line_ports().is_none()
            && !matches!(self.nodes[node].class(), Class::ExactLinear | Class::Smooth)
        {
            return Err(SystemError::Unsupported {
                node: self.names[node].clone(),
                needs: "a port tangent or the non-smooth solve tier",
            });
        }
        // Forward difference of the port map. The candidate slot is free to
        // overwrite: a component without a declared tangent is proposed
        // again, at the accepted inputs, before its block can be accepted.
        // The step is `sqrt(eps)` of the input's scale; where that moves an
        // output by less than its rounding, the column is differenced again
        // at a longer step (`resolve_column`).
        //
        // A perturbed proposal that asks for a shorter step is not the
        // step's: the component took the step at the iterate itself, and
        // the perturbation, the kernel's own, has carried one input over
        // the edge of an event (a slip that passes zero just beyond where
        // this step ends at rest). The column is then differenced on the
        // other side of the iterate, where the component still takes the
        // step. Only if it asks on both sides is the step shortened, to
        // the shorter of the two.
        let index = inputs.start + local;
        let base = self.inputs[index];
        let scale = base.abs().max(group_scale[self.input_group[index]]);
        let step = f64::EPSILON.sqrt() * if scale > 0.0 { scale } else { 1.0 };
        let mut asked: Option<f64> = None;
        for side in [step, -step] {
            self.inputs[index] = base + side;
            let proposal = self.propose_perturbed(node, dt_s, cx);
            self.inputs[index] = base;
            match proposal {
                Err(error) => return Err(self.refused(node, error, cx)),
                Ok(Proposal::Shorten { dt_s: shorter }) => {
                    let shorter = self.shorter(node, dt_s, shorter)?;
                    asked = match (asked, shorter) {
                        (Some(first), Some(second)) => Some(first.min(second)),
                        (first, second) => first.or(second),
                    };
                }
                Ok(Proposal::Step(_)) => {
                    for (column, current) in self.scratch_outputs[..width]
                        .iter_mut()
                        .zip(&self.outputs[outputs])
                    {
                        *column = (*column - current) / side;
                    }
                    // No step longer than `2^-10` of the input's scale, or
                    // of one unit where that is less.
                    let widest = scale.max(1.0) * WIDEST;
                    self.resolve_column(node, index, side, widest, dt_s, cx)?;
                    return Ok(None);
                }
            }
        }
        Ok(asked)
    }

    /// Propose `node` at the current inputs into `scratch_outputs`, on its
    /// load line if it is on one: a proposal at an input the kernel has
    /// moved, whose outputs difference a column of its port map.
    fn propose_perturbed(
        &mut self,
        node: usize,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Proposal, OwnerError> {
        let (inputs, width) = (self.node_inputs(node), self.node_outputs(node).len());
        let line = self.lines.get(self.line_of[node]).and_then(Option::as_ref);
        if let Some(line) = line {
            self.nodes[node].propose_on_line(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs],
                &line.compliance,
                self.candidate[node].as_mut(),
                &mut self.scratch_outputs[..width],
                cx,
            )
        } else {
            self.nodes[node].propose(
                self.accepted[node].as_ref(),
                dt_s,
                &self.inputs[inputs],
                self.candidate[node].as_mut(),
                &mut self.scratch_outputs[..width],
                cx,
            )
        }
    }

    /// Propose `node` with input `index` moved by `moved`, into
    /// `scratch_outputs`, and put the input back: whether the component
    /// takes the step there. One that asks to shorten or refuses there does
    /// not (the move is the kernel's own), unless the step is cancelled.
    fn steps_at(
        &mut self,
        node: usize,
        index: usize,
        moved: f64,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<bool, SystemError> {
        let base = self.inputs[index];
        self.inputs[index] = base + moved;
        let proposal = self.propose_perturbed(node, dt_s, cx);
        self.inputs[index] = base;
        match proposal {
            Ok(Proposal::Step(_)) => Ok(true),
            Err(_) if cx.is_cancel_requested() => Err(SystemError::Cancelled),
            Ok(Proposal::Shorten { .. }) | Err(_) => Ok(false),
        }
    }

    /// Resolve the column of `node` that `scratch_outputs` holds,
    /// differenced at `side` of input `index`, against the rounding of its
    /// outputs, the input moved no further than `widest`; the entries left
    /// unresolved are counted in the step's report.
    ///
    /// A step of `sqrt(eps)` of an input's scale can move an output by less
    /// than the output's rounding: the input is small where the output is
    /// not (a mass that moves with little force on it; a spring at rest
    /// that holds a load). The entry is then a count of rounding units, or
    /// zero: beside a 2 kg mass at 1 m/s under 1e-6 N a load line was
    /// handed no mobility at all. An entry is resolved when its change is
    /// at least [`RESOLVED`] of its output, or its output is zero. While one
    /// is not, the column is differenced again on the same side, at the step
    /// that aims the least resolved change at [`AIMED`] of its output (the
    /// change of a linear answer is in proportion to the step), at most
    /// [`LONGER_STEPS`] times. Each entry keeps the difference that first
    /// resolved it, so a column the first step resolves is the column it
    /// was, and an entry of a curved answer is not taken at a step longer
    /// than its own needs.
    ///
    /// An entry resolved at a longer step is differenced at the longest of
    /// them on the other side of the iterate, and kept if the two agree
    /// within [`AGREED`] of it: an output with a kink inside the longer step
    /// on either side gives two that do not, and the entry is left as the
    /// first step had it, never the secant across the kink. So is one that
    /// the component does not take on the other side (it asks to shorten
    /// there, or refuses). A longer step that changes an output by more than
    /// [`LOCAL`] of it has jumped, and the entry keeps what it had. An entry
    /// still unresolved at `widest`, or when the component does not take a
    /// longer step, keeps its latest difference. Each entry left unresolved
    /// is counted: an output that does not answer the input at all among
    /// them, its entry nothing.
    fn resolve_column(
        &mut self,
        node: usize,
        index: usize,
        side: f64,
        widest: f64,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<(), SystemError> {
        let outputs = self.node_outputs(node);
        let width = outputs.len();
        let mut longer = 1.0_f64;
        for ((entry, &column), &current) in self.scratch_entries[..width]
            .iter_mut()
            .zip(&self.scratch_outputs[..width])
            .zip(&self.outputs[outputs.clone()])
        {
            let change = (column * side).abs();
            *entry = if resolved(change, current) {
                Entry::First
            } else {
                longer = longer.max(aimed(change, current));
                Entry::Open
            };
        }
        if !self.scratch_entries[..width].contains(&Entry::Open) {
            return Ok(());
        }
        self.scratch_first[..width].copy_from_slice(&self.scratch_outputs[..width]);
        self.scratch_kept[..width].copy_from_slice(&self.scratch_outputs[..width]);
        let mut step = side.abs();
        let mut moved = 0.0_f64;
        for _ in 0..LONGER_STEPS {
            if step >= widest || !self.scratch_entries[..width].contains(&Entry::Open) {
                break;
            }
            step = (step * longer).min(widest);
            if !self.steps_at(node, index, step.copysign(side), dt_s, cx)? {
                break;
            }
            moved = step.copysign(side);
            longer = 1.0;
            for (((kept, entry), &there), &current) in self.scratch_kept[..width]
                .iter_mut()
                .zip(&mut self.scratch_entries[..width])
                .zip(&self.scratch_outputs[..width])
                .zip(&self.outputs[outputs.clone()])
            {
                if *entry != Entry::Open {
                    continue;
                }
                let change = there - current;
                if change.abs() > LOCAL * current.abs() {
                    *entry = Entry::Left;
                    continue;
                }
                *kept = change / moved;
                if resolved(change.abs(), current) {
                    *entry = Entry::Taken;
                } else {
                    longer = longer.max(aimed(change.abs(), current));
                }
            }
        }
        if self.scratch_entries[..width].contains(&Entry::Taken) {
            let stepped = self.steps_at(node, index, -moved, dt_s, cx)?;
            for ((((kept, entry), &first), &there), &current) in self.scratch_kept[..width]
                .iter_mut()
                .zip(&mut self.scratch_entries[..width])
                .zip(&self.scratch_first[..width])
                .zip(&self.scratch_outputs[..width])
                .zip(&self.outputs[outputs.clone()])
            {
                if *entry == Entry::Taken
                    && !(stepped
                        && ((there - current) / -moved - *kept).abs() <= AGREED * kept.abs())
                {
                    *kept = first;
                    *entry = Entry::Left;
                }
            }
        }
        let left = self.scratch_entries[..width]
            .iter()
            .filter(|&&entry| matches!(entry, Entry::Open | Entry::Left))
            .count();
        self.unresolved = self
            .unresolved
            .saturating_add(u32::try_from(left).unwrap_or(u32::MAX));
        self.scratch_outputs[..width].copy_from_slice(&self.scratch_kept[..width]);
        Ok(())
    }

    /// The largest power [W] crossing a power port of `nodes` at the current
    /// inputs and outputs: the absolute term sum of its power pairing.
    fn block_power_scale(&self, nodes: &[usize]) -> f64 {
        let mut largest = 0.0_f64;
        for &node in nodes {
            let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
            for (port, spec) in self.nodes[node].ports().iter().enumerate() {
                if !matches!(spec.carrier(), Carrier::Power(_)) {
                    continue;
                }
                let received = self.layouts[node].input_range(port);
                let returned = self.layouts[node].output_range(port);
                let power = power_pairing_scale(
                    spec,
                    &self.inputs[inputs.start + received.start..inputs.start + received.end],
                    &self.outputs[outputs.start + returned.start..outputs.start + returned.end],
                );
                largest = largest.max(power);
            }
        }
        largest
    }

    /// The second convergence test of a block: the residuals carried once
    /// through `slopes`. `None` when no wire's target moves by more than
    /// its own tolerance,
    /// `sum_j |t_ij| max(0, |r_j| - f_j) <= rel S_i + f_i` for every wire
    /// `i`; otherwise the open input (a column) that contributes most,
    /// over all its entries, to the wire that is furthest beyond it.
    ///
    /// `f` is a wire's rounding floor at the current inputs: the absolute
    /// tolerance plus `ROUNDING_UNITS` of `eps` times the larger of the
    /// wire's own value and the terms its driver sums, plus the size of
    /// the inputs its target depends on through the slopes (for an energy
    /// rate also the block's largest power), plus `kind_rounding`, which
    /// is zero outside a block with a load line. `S` is `scale`. `floor`
    /// and `carried` are scratch of the block's size.
    #[allow(clippy::too_many_arguments)]
    fn carried_beyond(
        &self,
        unknowns: &[usize],
        residual: &[f64],
        slopes: &[(usize, usize, f64)],
        scale: &[f64],
        kind_rounding: &[f64],
        power_scale: f64,
        floor: &mut [f64],
        carried: &mut [f64],
    ) -> Option<usize> {
        floor.fill(0.0);
        carried.fill(0.0);
        for &(row, column, through) in slopes {
            floor[row] += (through * self.inputs[unknowns[column]]).abs();
        }
        for (row, &index) in unknowns.iter().enumerate() {
            let own = self.inputs[index]
                .abs()
                .max(self.drivers[index].magnitude(&self.outputs, &self.boundary_values));
            floor[row] =
                self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * (floor[row] + own);
            if self.input_is_energy_rate[index] {
                floor[row] += ROUNDING_UNITS * f64::EPSILON * power_scale;
            }
            floor[row] += kind_rounding[row];
        }
        for &(row, column, through) in slopes {
            let above = residual[column].abs() - floor[column];
            if above > 0.0 {
                carried[row] += through.abs() * above;
            }
        }
        // The wire whose target is moved furthest beyond its tolerance. A
        // value that is not a number fails the comparison and is taken.
        let mut moved = None;
        let mut furthest = 0.0_f64;
        for row in 0..unknowns.len() {
            let allowed = self.config.relative_tolerance * scale[row] + floor[row];
            let excess = carried[row] - allowed;
            if !(excess <= 0.0) && !(excess <= furthest) {
                moved = Some(row);
                furthest = excess;
            }
        }
        let moved = moved?;
        // The open input that does most of it, all its entries on the
        // moved wire together; the moved wire itself if nothing can be
        // ranked. `carried` is free again and is now per input.
        carried.fill(0.0);
        for &(row, column, through) in slopes {
            if row == moved {
                carried[column] += through.abs() * (residual[column].abs() - floor[column]);
            }
        }
        let mut open = moved;
        let mut largest = 0.0_f64;
        for (column, &part) in carried.iter().enumerate() {
            if part > largest {
                largest = part;
                open = column;
            }
        }
        Some(open)
    }

    /// The refusal of a block that did not converge, naming input `row` of
    /// its unknowns.
    fn not_converged(
        &self,
        nodes: &[usize],
        unknowns: &[usize],
        residual: &[f64],
        iterations: u32,
        row: usize,
    ) -> SystemError {
        SystemError::NotConverged {
            nodes: self.names_of(nodes),
            iterations,
            worst: self.input_label[unknowns[row]].clone(),
            residual: residual[row],
        }
    }

    /// Solve one looped block by Newton iteration on its interface inputs.
    fn solve_block(
        &mut self,
        block: usize,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        // The block's lists and workspace are lent to the solve and handed
        // back, so a step allocates nothing for them.
        let unknowns = core::mem::take(&mut self.blocks[block].unknowns);
        let nodes = core::mem::take(&mut self.blocks[block].nodes);
        let mut work = core::mem::take(&mut self.block_work[block]);
        let result = self.solve_lent(block, &unknowns, &nodes, &mut work, dt_s, cx);
        self.blocks[block].unknowns = unknowns;
        self.blocks[block].nodes = nodes;
        self.block_work[block] = work;
        if matches!(result, Ok(Ok(_))) {
            // The flow inputs of a load line held released flows during
            // the solve; what is published is what crossed.
            self.line_flows(block);
        }
        result
    }

    /// [`Self::solve_block`] on a block's own lists and workspace.
    #[allow(clippy::too_many_lines)]
    fn solve_lent(
        &mut self,
        block: usize,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        for &index in unknowns {
            let start = self.accepted_inputs[index];
            if self.input_is_temperature[index] && !(start.is_finite() && start > 0.0) {
                return Err(SystemError::GuessRequired {
                    port: self.input_label[index].clone(),
                });
            }
            self.inputs[index] = start;
        }
        work.group_scale.resize(self.group_count, 0.0);
        work.rounding.fill(0.0);
        work.allowance.fill(0.0);
        // A block with a component on a load line: the line of this step,
        // found from the starting values, before anything is judged.
        if self.lines[block].is_some()
            && let Some(shorter) = self.find_line(block, unknowns, nodes, work, dt_s, cx)?
        {
            return Ok(Err(shorter));
        }
        let BlockWork {
            residual,
            rounding,
            allowance,
            group_scale,
            gradient,
            condensed,
            column_readers,
            readers,
            slots,
            slot_start,
            varying,
            varying_terms,
            affine,
            kept_dt_s,
            slopes,
            floor,
            carried,
            scale,
            kind_rounding,
            found,
        } = work;
        let affine = *affine;
        // Until the first update: the inputs of every component other than
        // the line component are those finding the line proposed it at,
        // unless shown otherwise below (`Found::repeats`).
        let mut first_pass = true;
        // Whether this step has already used the factors an affine block
        // kept from an earlier one.
        let mut kept_used = false;
        // Whether this step has assembled the block's Jacobian yet, so that
        // `slopes` are this step's own.
        let mut assembled = false;
        // A trial starts from the values the last step ended with. For the
        // component on a load line those are the flows of the last step, not
        // the released flows of this one, and whatever it answers there is
        // answered from the wrong flows: a slider that enters a step with a
        // slip its ports were not told of would see a stop inside every
        // step, and the shorter step it asks for would start from the same
        // values and ask again, without end. So before it is first asked in
        // a trial, its inputs (the line's flows and ends, and any other,
        // such as a span) take one update with its efforts held, which is
        // exact where the other components are affine; every other unknown
        // keeps its start.
        if self.lines[block].is_some() {
            self.release(block, unknowns, residual, condensed)?;
        }
        let mut iterations = 0_u32;
        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            for (position, &node) in nodes.iter().enumerate() {
                // In the first pass, a component that finding the line
                // proposed, at bit for bit the inputs it holds now, is not
                // proposed again: its accepted state and the step are the
                // same (states change only when a step is published), and a
                // proposal is pure (fs-ports: equal arguments give
                // bitwise-equal results), so its candidate, outputs and
                // receipt are already this proposal's. Only a component
                // that declares its tangent: differencing a component
                // overwrites its candidate. The line component's inputs
                // were moved by `release`; it was not proposed there.
                if first_pass
                    && !found.proposed.is_empty()
                    && found.repeats(position, &self.inputs[self.node_inputs(node)])
                {
                    continue;
                }
                if let Some(shorter) = self.propose(node, dt_s, cx)? {
                    return Ok(Err(shorter));
                }
            }
            // The scale of an unknown is the size of the terms its driver
            // sums, not of their total: a force balance near rest is a small
            // difference of large efforts, and its residual cannot be driven
            // below the roundoff of those efforts.
            //
            // A flow on a load line is the exception. Its input holds the
            // released flow and its driver sums what the line component's
            // own efforts take away, so by those sizes a slider that holds
            // a mass at rest would be judged against the flow it would
            // have if it let go. Its size is what crosses and what the
            // other components give.
            group_scale.fill(0.0);
            for (row, &index) in unknowns.iter().enumerate() {
                let driver = &self.drivers[index];
                let target = driver.value(&self.outputs, &self.boundary_values);
                residual[row] = self.inputs[index] - target;
                let group = self.input_group[index];
                let (own, terms) = match self.crossing(block, row) {
                    Some((flow, _)) => (
                        flow.abs(),
                        driver.magnitude_of_others(&self.outputs, &self.boundary_values),
                    ),
                    None => (
                        self.inputs[index].abs(),
                        driver.magnitude(&self.outputs, &self.boundary_values),
                    ),
                };
                group_scale[group] = group_scale[group].max(own).max(terms);
                // The wire's own size: what a block with a load line
                // judges it by (below).
                scale[row] = own.max(terms);
            }
            // A residual cannot be driven below the rounding of the values
            // it is computed from. When a component's output is a small
            // difference of large inputs (a heat rate between two bodies near
            // one temperature, a spring force between two ends that move
            // together), that rounding is set by the inputs, not by the
            // output, and no tolerance relative to the output covers it.
            //
            // That estimate is first order. A heat rate that a component
            // forms as a product of an effort and a flow (a damper's force
            // times its rate) is rounded at second order when both factors
            // are themselves rounding: a body held at rest by forces that
            // cancel moves at the rounding of those forces, and its damper
            // then reports the square of that. No first-order estimate sees
            // it, and the heat has no size of its own to be relative to. The
            // largest power crossing a power port of the block bounds it (it
            // has the same rounding-sized flow against an effort that is not
            // rounding), so an energy rate is also allowed rounding in units
            // of that power.
            let power_scale = self.block_power_scale(nodes);
            let mut worst = (0_usize, 0.0_f64);
            let mut converged = true;
            let on_a_line = self.lines[block].is_some();
            for (row, &index) in unknowns.iter().enumerate() {
                // The interface audit is told only this part: it keeps its
                // own relative tolerance, so a loose solver tolerance is
                // still caught there.
                allowance[row] =
                    self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * rounding[row];
                if self.input_is_energy_rate[index] {
                    allowance[row] += ROUNDING_UNITS * f64::EPSILON * power_scale;
                }
                // In a block with a load line every wire is relative to
                // its own size, not to the largest of its kind in the
                // block. The line component is differenced, so an update
                // leaves each wire with what the differencing left of the
                // step's change, and the audit of a junction knows that
                // junction only and is told rounding, never the relative
                // tolerance here. Taken by the block's scale, a mass held
                // beside one that moves was refused, and a block riding a
                // cart when the friction between them passed through
                // zero.
                //
                // A wire's own size is not its floor, though. The block's
                // unknowns are solved together, and the rounding of a
                // large wire is in the small ones: the last digit of a
                // force balance of 11.7 N moves a wire of 1e-25 N beside
                // it by a thousandth of itself from one iterate to the
                // next. So a wire that is not a flow on the line is
                // allowed the rounding of the largest wire of its kind in
                // the block, and the audit is told that for it.
                //
                // A flow on the line has the rounding of the difference it
                // is, which no assembly estimates, and the solve takes it
                // only at half of what the audit is then told, so that a
                // step taken at the edge of the one is not refused at the
                // edge of the other.
                kind_rounding[row] = 0.0;
                let mut taken_at = allowance[row];
                if let Some((_, rounding)) = self.crossing(block, row) {
                    allowance[row] += rounding;
                    taken_at = 0.5 * allowance[row];
                } else if on_a_line {
                    kind_rounding[row] =
                        ROUNDING_UNITS * f64::EPSILON * group_scale[self.input_group[index]];
                    allowance[row] += kind_rounding[row];
                    taken_at = allowance[row];
                } else {
                    scale[row] = group_scale[self.input_group[index]];
                }
                let allowed = self.config.relative_tolerance * scale[row] + taken_at;
                let excess = residual[row].abs() - allowed;
                if !(excess <= 0.0) {
                    converged = false;
                    if !(excess <= worst.1) {
                        worst = (row, excess);
                    }
                }
            }
            // A residual that passes its own test can still be large for
            // what reads it. An effort is a gauge quantity: a temperature
            // open by 1e-12 of 300 K is 3e-10 K, and a conductance that
            // reads it returns a heat rate that is wrong by that much times
            // its conductance, while the heat wire agrees with what the
            // conductance returned and passes too. So the residuals are
            // also carried once through the components that read them, by
            // the slopes of the block's Jacobian (`carried_beyond`), and
            // may not move any wire's target by more than that wire's own
            // tolerance.
            //
            // Only the part of a residual above its own rounding floor is
            // carried. A force balance near rest sits at the rounding of
            // forces that cancel, which is far above the rounding of its
            // small net value and which nothing can lower; carried whole,
            // it would keep a body's velocity wire from ever closing. The
            // floor of a wire is therefore the rounding of three things:
            // the terms its driver sums, its own value, and the inputs its
            // target depends on through the slopes. (The first test has
            // the first two through its scale; the allowance it adds has
            // only the third.)
            //
            // What this bounds for each wire is its residual and the
            // first-order image of the residuals under one more pass
            // through the block. It does not bound the whole first-order
            // error `(I - T)^-1 r` of a loop whose gain exceeds one.
            //
            // The slopes are always this step's own. What a block keeps
            // from an earlier step may not decide a step: a checkpoint does
            // not carry it and a discarded trial does not roll it back, so
            // a verdict that read it would differ between a run and its
            // replay. At a step's first check there are none yet, and they
            // are formed below at the iterate being judged; at a later
            // check they are those of the assembly that produced the
            // update, one iterate earlier.
            let mut judged_below = false;
            if converged {
                if residual.iter().all(|value| *value == 0.0) {
                    // An exact fixed point needs no slopes.
                } else if assembled {
                    if let Some(open) = self.carried_beyond(
                        unknowns,
                        residual,
                        slopes,
                        scale,
                        kind_rounding,
                        power_scale,
                        floor,
                        carried,
                    ) {
                        converged = false;
                        worst.0 = open;
                    }
                } else {
                    converged = false;
                    judged_below = true;
                }
            }
            if converged {
                for (row, &index) in unknowns.iter().enumerate() {
                    self.input_allowance[index] = allowance[row];
                }
                return Ok(Ok(iterations));
            }
            if iterations >= self.config.max_iterations && !judged_below {
                return Err(self.not_converged(nodes, unknowns, residual, iterations, worst.0));
            }

            // J = I - (d driver / d outputs) (d outputs / d inputs).
            rounding.fill(0.0);
            if affine && !kept_used && kept_dt_s.to_bits() == dt_s.to_bits() {
                // Every part of the block is affine in its inputs, so the
                // Jacobian is the one factored at this step length before.
                // Only the rounding estimate depends on the inputs. If a
                // part's slopes moved with its state after all, the update
                // does not converge and the next iteration derives them.
                kept_used = true;
                for &(row, column, through) in slopes.iter() {
                    rounding[row] += (through * self.inputs[unknowns[column]]).abs();
                }
            } else {
                *kept_dt_s = f64::NAN;
                slopes.clear();
                condensed.clear();
                // The slopes of the drivers that are not linear, at the
                // outputs of this iterate.
                for (row, terms) in varying.iter() {
                    self.drivers[unknowns[*row]].gradient(&self.outputs, gradient);
                    for (&reader, &(_, slope)) in
                        varying_terms[terms.clone()].iter().zip(gradient.iter())
                    {
                        if reader != usize::MAX {
                            readers[reader].2 = slope;
                        }
                    }
                }
                for (column, &index) in unknowns.iter().enumerate() {
                    let node = self.input_owner[index];
                    let local = index - self.input_start[node];
                    // The first assembly takes a repeated component's
                    // column from finding the line, at the same inputs.
                    if first_pass
                        && !found.proposed.is_empty()
                        && let Some(kept) = found.column_of(column, self.node_outputs(node).len())
                    {
                        self.scratch_outputs[..kept.len()].copy_from_slice(kept);
                    } else if let Some(shorter) =
                        self.tangent_column(node, local, dt_s, group_scale, cx)?
                    {
                        return Ok(Err(shorter));
                    }
                    // Every row that reads this component's outputs, where
                    // the plan keeps its entry in this column.
                    let reading = &readers[column_readers[column].clone()];
                    for (&(row, offset, slope), &slot) in
                        reading.iter().zip(&slots[slot_start[column]..])
                    {
                        let through = slope * self.scratch_outputs[offset];
                        condensed.subtract(slot, through);
                        rounding[row] += (through * self.inputs[index]).abs();
                        slopes.push((row, column, through));
                    }
                }
                if !condensed.factor() {
                    return Err(SystemError::Singular {
                        nodes: self.names_of(nodes),
                    });
                }
                if affine {
                    *kept_dt_s = dt_s;
                }
            }
            assembled = true;
            if judged_below {
                // The step's first check, judged with the slopes just
                // formed at its own iterate.
                let Some(open) = self.carried_beyond(
                    unknowns,
                    residual,
                    slopes,
                    scale,
                    kind_rounding,
                    power_scale,
                    floor,
                    carried,
                ) else {
                    // Taken with no update. Differencing a component wrote
                    // a perturbed proposal to its candidate slot, so each
                    // component without a declared tangent is proposed
                    // once more at the accepted inputs. (A declared
                    // tangent writes to scratch only, and kept factors
                    // call no component.)
                    for &node in nodes {
                        if self.nodes[node].has_tangent() {
                            continue;
                        }
                        if let Some(shorter) = self.propose(node, dt_s, cx)? {
                            return Ok(Err(shorter));
                        }
                    }
                    // The audit is told the rounding of this iterate, which
                    // the assembly above has just estimated. The allowance
                    // the first test used was formed before it and has
                    // none; told that, the audit of a junction whose flow
                    // has gone to zero (a body come to the temperature of
                    // what it touches) would allow a rounding-sized
                    // difference of large terms nothing at all.
                    for (row, &index) in unknowns.iter().enumerate() {
                        let mut told = self.config.absolute_tolerance
                            + ROUNDING_UNITS * f64::EPSILON * rounding[row];
                        if self.input_is_energy_rate[index] {
                            told += ROUNDING_UNITS * f64::EPSILON * power_scale;
                        }
                        if let Some((_, rounding)) = self.crossing(block, row) {
                            told += rounding;
                        }
                        told += kind_rounding[row];
                        self.input_allowance[index] = told;
                    }
                    return Ok(Ok(iterations));
                };
                if iterations >= self.config.max_iterations {
                    return Err(self.not_converged(nodes, unknowns, residual, iterations, open));
                }
            }
            iterations += 1;
            first_pass = false;
            for value in residual.iter_mut() {
                *value = -*value;
            }
            if !condensed.solve(residual) {
                return Err(SystemError::Singular {
                    nodes: self.names_of(nodes),
                });
            }
            for (row, &index) in unknowns.iter().enumerate() {
                self.inputs[index] += residual[row];
            }
        }
    }

    /// Fill every candidate for a step of `dt_s`.
    fn trial(&mut self, dt_s: f64, cx: &Cx<'_>) -> Result<Trial, SystemError> {
        self.input_allowance.fill(0.0);
        let t0 = self.time.value();
        for ((_, value), slot) in self.boundaries.iter().zip(&mut self.boundary_values) {
            value(t0, dt_s, slot);
        }
        let mut iterations = 0_u32;
        for block in 0..self.blocks.len() {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            // Inputs that come from earlier blocks, boundaries and constants.
            for position in 0..self.blocks[block].known.len() {
                let index = self.blocks[block].known[position];
                self.inputs[index] =
                    self.drivers[index].value(&self.outputs, &self.boundary_values);
            }
            if self.blocks[block].looped() {
                match self.solve_block(block, dt_s, cx)? {
                    Ok(spent) => iterations += spent,
                    Err(shorter) => return Ok(Trial::Shorten { dt_s: shorter }),
                }
            } else {
                for position in 0..self.blocks[block].nodes.len() {
                    let node = self.blocks[block].nodes[position];
                    if let Some(shorter) = self.propose(node, dt_s, cx)? {
                        return Ok(Trial::Shorten { dt_s: shorter });
                    }
                }
            }
        }
        Ok(Trial::Done { iterations })
    }

    /// Audit every candidate and fill the per-port energies.
    // `node` indexes six parallel per-component tables, not only `defects`.
    #[allow(clippy::needless_range_loop)]
    fn audit(
        &mut self,
        dt_s: f64,
        defects: &mut [f64],
        mass_audits: &mut [Option<MassStepAudit>],
    ) -> Result<(), SystemError> {
        for node in 0..self.nodes.len() {
            let measure = |state: &dyn Any| {
                self.nodes[node]
                    .stored_energy_j(state)
                    .map_err(|error| SystemError::Component {
                        node: self.names[node].clone(),
                        phase: Phase::Measure,
                        error,
                    })
            };
            let stored_old = measure(self.accepted[node].as_ref())?;
            let stored_new = measure(self.candidate[node].as_ref())?;
            let (inputs, outputs) = (self.node_inputs(node), self.node_outputs(node));
            let ports = self.nodes[node].ports();
            let audit = audit_step(
                ports,
                &self.layouts[node],
                dt_s,
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
            defects[node] = audit.balance_defect_j;
            for (port, spec) in ports.iter().enumerate() {
                let input =
                    &self.inputs[inputs.start..inputs.end][self.layouts[node].input_range(port)];
                let output = &self.outputs[outputs.start..outputs.end]
                    [self.layouts[node].output_range(port)];
                let global = self.port_start[node] + port;
                self.port_energy[global] = spec.energy_in_j(dt_s, input, output);
                self.port_entropy[global] = spec.entropy_in_j_per_k(dt_s, input, output);
                self.port_mass[global] = spec.mass_in_kg(dt_s, input, output);
                for (carried, value) in self.port_carried[global].iter_mut().enumerate() {
                    *value = spec.carried_quantity_in(dt_s, input, output, carried);
                }
            }
            if self.ledger.nodes[node].has_streams {
                let measure_mass = |state: &dyn Any| {
                    self.nodes[node]
                        .stored_mass_kg(state)
                        .map_err(|error| SystemError::Component {
                            node: self.names[node].clone(),
                            phase: Phase::Measure,
                            error,
                        })
                };
                let old_mass = measure_mass(self.accepted[node].as_ref())?;
                let new_mass = measure_mass(self.candidate[node].as_ref())?;
                mass_audits[node] = audit_mass_step(
                    ports,
                    &self.layouts[node],
                    dt_s,
                    &self.inputs[inputs.clone()],
                    &self.outputs[outputs.clone()],
                    old_mass,
                    new_mass,
                )
                .map_err(|violation| SystemError::Contract {
                    node: self.names[node].clone(),
                    violation,
                })?;
            }
            for carried in 0..self.carried_audits[node].len() {
                let measure = |state: &dyn Any| self.nodes[node]
                    .stored_carried_quantity(state, carried)
                    .map_err(|error| SystemError::Component {
                        node: self.names[node].clone(), phase: Phase::Measure, error,
                    });
                let old = measure(self.accepted[node].as_ref())?;
                let new = measure(self.candidate[node].as_ref())?;
                self.carried_audits[node][carried] = audit_carried_step(
                    ports, &self.layouts[node], dt_s,
                    &self.inputs[inputs.clone()], &self.outputs[outputs.clone()], carried, old, new,
                ).map_err(|violation| SystemError::Contract {
                    node: self.names[node].clone(), violation,
                })?;
            }
        }
        Ok(())
    }

    /// Advance by at most `dt_max_s`.
    ///
    /// The step is the largest that every component's `limit_dt_s` and
    /// `next_event_s` allow, shortened further if a component asks during the
    /// trial. All candidates are audited and admitted before any is
    /// published.
    ///
    /// # Errors
    /// Any [`SystemError`]; the system is unchanged when a step fails.
    #[allow(clippy::too_many_lines)]
    pub fn step(&mut self, dt_max_s: f64, cx: &Cx<'_>) -> Result<StepReport, SystemError> {
        if !(dt_max_s.is_finite() && dt_max_s > 0.0) {
            return Err(SystemError::Step { dt_s: dt_max_s });
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
        let mut dt_s = dt_max_s;
        for node in 0..self.nodes.len() {
            let limited = self.nodes[node]
                .limit_dt_s(self.accepted[node].as_ref(), dt_s)
                .map_err(|error| SystemError::Component {
                    node: self.names[node].clone(),
                    phase: Phase::Negotiate,
                    error,
                })?;
            // The limit and the event are checked one at a time: `min` would
            // silently drop a NaN from either.
            let refuse = |returned_s: f64| SystemError::Negotiation {
                node: self.names[node].clone(),
                proposed_s: dt_s,
                returned_s,
            };
            if !(limited > 0.0 && limited <= dt_s) {
                return Err(refuse(limited));
            }
            let mut next = limited;
            if let Some(event) = self.nodes[node].next_event_s(self.accepted[node].as_ref()) {
                if !(event > 0.0) {
                    return Err(refuse(event));
                }
                next = next.min(event);
            }
            dt_s = next;
        }

        self.unresolved = 0;
        let mut shortenings = 0_u32;
        let iterations = loop {
            match self.trial(dt_s, cx)? {
                Trial::Done { iterations } => break iterations,
                // `shorter` was checked against `dt_s` where it was asked for.
                Trial::Shorten { dt_s: shorter } => {
                    shortenings += 1;
                    if shortenings > self.config.max_shortenings {
                        return Err(SystemError::StepUnderflow { dt_s: shorter });
                    }
                    dt_s = shorter;
                }
            }
        };

        // The step's lists are taken from the system and given back when
        // the step is published, so an accepted step allocates nothing for
        // them. A refused step drops them and the next one makes new ones.
        let StepScratch {
            mut defects,
            mut mass_audits,
            mut members,
            mut group_members,
            mut connection_rows,
            mut sink_heat,
        } = core::mem::take(&mut self.step_scratch);
        defects.clear();
        defects.resize(self.nodes.len(), 0.0);
        mass_audits.clear();
        mass_audits.resize(self.nodes.len(), None);
        self.audit(dt_s, &mut defects, &mut mass_audits)?;

        // Connection and boundary rows of this step.
        members.clear();
        connection_rows.clear();
        sink_heat.clear();
        sink_heat.resize(self.ledger.sinks.len(), 0.0);
        let (mut generated, mut entropy_scale) = (0.0_f64, 0.0_f64);
        let mut entropy_slack = 0.0_f64;
        for receipt in &self.receipts {
            generated += receipt.entropy_generated_j_per_k;
            entropy_scale += receipt.entropy_generated_j_per_k;
        }
        for (index, group) in self.audit_groups.iter().enumerate() {
            // A heat route that shares a receiver with an earlier one is
            // audited in that one's row.
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
            let mass_scale = members
                .iter()
                .map(|port| self.port_mass[self.port_start[port.node] + port.port].abs())
                .fold(0.0_f64, f64::max);
            let (mut energy, mut magnitude, mut entropy) = (0.0_f64, 0.0_f64, 0.0_f64);
            let mut mass = Sum::default();
            let mut mass_magnitude = 0.0_f64;
            for port in &members {
                let global = self.port_start[port.node] + port.port;
                energy += self.port_energy[global];
                if mass_scale > 0.0 {
                    let scaled_mass = self.port_mass[global] / mass_scale;
                    mass.add(scaled_mass);
                    mass_magnitude += scaled_mass.abs();
                }
                let spec = &self.nodes[port.node].ports()[port.port];
                magnitude += if spec.field_metric().is_some() {
                    dt_s * power_pairing_scale(
                        spec,
                        &self.inputs[self.input_range(*port)],
                        &self.outputs[self.output_range(*port)],
                    )
                } else {
                    self.port_energy[global].abs()
                };
                entropy += self.port_entropy[global];
                entropy_scale += self.port_entropy[global].abs();
            }
            for &connection in group {
                if let Connection::HeatRoute {
                    source,
                    remainder: Some((sink, fraction, kelvin)),
                    ..
                } = &self.connections[connection]
                {
                    let emitted = -self.port_energy[self.port_start[source.node] + source.port];
                    let heat = fraction * emitted;
                    sink_heat[*sink] += heat;
                    energy += heat;
                    magnitude += heat.abs();
                    entropy += heat / kelvin;
                    entropy_scale += (heat / kelvin).abs();
                }
            }
            // What rounding and the configured absolute tolerance let the
            // block solve leave in each input, as energy: a power port's
            // input pairs with its output, a heat receiver's input is the
            // heat rate itself, and of the rates a storage receives from a
            // stream the enthalpy rate is the energy.
            let mut solve_slack_j = 0.0_f64;
            let mut mass_slack_kg = 0.0_f64;
            for port in &members {
                let spec = &self.nodes[port.node].ports()[port.port];
                let outputs = self.output_range(*port);
                for (offset, index) in self.input_range(*port).enumerate() {
                    if spec.carrier() == Carrier::Stream
                        && spec.role() == PortRole::FlowIn
                        && offset == STREAM_MASS_RATE
                    {
                        mass_slack_kg += dt_s * self.input_allowance[index];
                    }
                    let slack_j = match (spec.carrier(), spec.role()) {
                        (Carrier::Heat, PortRole::FlowIn) => dt_s * self.input_allowance[index],
                        (Carrier::Stream, PortRole::FlowIn) if offset == STREAM_ENTHALPY_RATE => {
                            dt_s * self.input_allowance[index]
                        }
                        (Carrier::Power(_), _) => power_input_slack_j(
                            spec,
                            &self.outputs[outputs.clone()],
                            offset,
                            dt_s,
                            self.input_allowance[index],
                        ),
                        _ => 0.0,
                    };
                    solve_slack_j += slack_j;
                    // The same heat, left uncancelled at the temperature
                    // its receiver returned, is entropy the solve was
                    // allowed to leave unaccounted.
                    if spec.carrier() == Carrier::Heat && spec.role() == PortRole::FlowIn {
                        entropy_slack += slack_j / self.outputs[outputs.start + offset];
                    }
                }
            }
            let tolerance_j = self.config.interface_energy_tolerance * magnitude
                + solve_slack_j
                + f64::MIN_POSITIVE;
            if !(energy.is_finite() && tolerance_j.is_finite() && energy.abs() <= tolerance_j) {
                return Err(SystemError::Interface {
                    connection: index,
                    defect_j: energy,
                    tolerance_j,
                });
            }
            generated += entropy;
            let mass_defect_kg = mass.value() * mass_scale;
            let mass_tolerance_kg = ((128.0 * f64::EPSILON + self.config.relative_tolerance)
                * mass_magnitude)
                * mass_scale
                + mass_slack_kg;
            if !(mass_defect_kg.is_finite()
                && mass_tolerance_kg.is_finite()
                && mass_defect_kg.abs() <= mass_tolerance_kg)
            {
                return Err(SystemError::InterfaceMass {
                    connection: index,
                    defect_kg: mass_defect_kg,
                    tolerance_kg: mass_tolerance_kg,
                });
            }
            connection_rows.push((energy, entropy, mass_defect_kg));
            for carried in 0..self.connection_carried[index].len() {
                let amount = |port: &PortRef| self.port_carried[self.port_start[port.node] + port.port]
                    .get(carried).copied().unwrap_or(0.0);
                let scale = members.iter().map(|port| amount(port).abs()).fold(0.0_f64, f64::max);
                let mut sum = Sum::default();
                let mut magnitude = 0.0;
                let mut slack = 0.0;
                for port in &members {
                    if scale > 0.0 {
                        let value = amount(port) / scale;
                        sum.add(value);
                        magnitude += value.abs();
                    }
                    let spec = &self.nodes[port.node].ports()[port.port];
                    if spec.role() == PortRole::FlowIn && carried < spec.carried_quantities() {
                        slack += dt_s * self.input_allowance[
                            self.input_range(*port).start + STREAM_CARRIED + carried
                        ];
                    }
                }
                let defect = sum.value() * scale;
                let tolerance = (128.0 * f64::EPSILON + self.config.relative_tolerance)
                    * magnitude * scale + slack;
                if !(defect.is_finite() && tolerance.is_finite() && defect.abs() <= tolerance) {
                    return Err(SystemError::InterfaceCarried {
                        connection: index, carried, defect, tolerance,
                    });
                }
                self.connection_carried[index][carried] = defect;
            }
        }
        if !(generated >= -1.0e-9 * entropy_scale - entropy_slack) {
            return Err(SystemError::SecondLaw {
                generated_j_per_k: generated,
            });
        }

        for node in 0..self.nodes.len() {
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
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;

        // Publish: every state at once, then the clock and the ledger.
        for (accepted, candidate) in self.accepted.iter_mut().zip(&mut self.candidate) {
            core::mem::swap(accepted, candidate);
        }
        self.accepted_inputs.copy_from_slice(&self.inputs);
        self.accepted_outputs.copy_from_slice(&self.outputs);
        self.accepted_allowance.copy_from_slice(&self.input_allowance);
        self.time.add(dt_s);
        self.steps += 1;
        self.ledger.steps += 1;
        self.ledger.duration_s.add(dt_s);
        for (node, account) in self.ledger.nodes.iter_mut().enumerate() {
            let receipt = &self.receipts[node];
            let ports = self.port_start[node]..self.port_start[node + 1];
            account
                .port_energy_in_j
                .add(self.port_energy[ports].iter().sum());
            account.stored_change_j.add(receipt.stored_change_j);
            account.dissipated_j.add(receipt.dissipated_j);
            account.numerical_loss_j.add(receipt.numerical_loss_j);
            account
                .entropy_generated_j_per_k
                .add(receipt.entropy_generated_j_per_k);
            account.balance_defect_j.add(defects[node]);
            if let Some(audit) = &mass_audits[node] {
                let mass = account.mass.get_or_insert_with(Default::default);
                mass.port_mass_in_kg.add(audit.port_mass_in_kg);
                mass.stored_change_kg.add(audit.stored_change_kg);
                mass.balance_defect_kg.add(audit.balance_defect_kg);
                mass.audited_steps += 1;
            }
            for (account, audit) in account.carried.iter_mut().zip(&self.carried_audits[node]) {
                if let Some(audit) = audit {
                    let account = account.get_or_insert_with(Default::default);
                    account.port_in.add(audit.port_in);
                    account.stored_change.add(audit.stored_change);
                    account.balance_defect.add(audit.balance_defect);
                    account.audited_steps += 1;
                }
            }
        }
        for (account, &(energy, entropy, mass)) in
            self.ledger.connections.iter_mut().zip(&connection_rows)
        {
            account.energy_defect_j.add(energy);
            account.entropy_generated_j_per_k.add(entropy);
            account.mass_defect_kg.add(mass);
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
        self.step_scratch = StepScratch {
            defects,
            mass_audits,
            members,
            group_members,
            connection_rows,
            sink_heat,
        };
        for (account, (port, _)) in self.ledger.boundaries.iter_mut().zip(&self.boundaries) {
            let global = self.port_start[port.node] + port.port;
            account.energy_in_j.add(self.port_energy[global]);
            account.entropy_in_j_per_k.add(self.port_entropy[global]);
            account.mass_in_kg.add(self.port_mass[global]);
            for (total, &value) in account.carried_in.iter_mut().zip(&self.port_carried[global]) {
                total.add(value);
            }
        }
        Ok(StepReport {
            dt_s,
            time_s: self.time.value(),
            iterations,
            shortenings,
            unresolved: self.unresolved,
        })
    }
}

/// Convenience: a [`HeatSink`] from a name and a temperature.
///
/// # Errors
/// An invalid name.
pub fn heat_sink(name: &str, temperature_k: f64) -> Result<HeatSink, fs_ports::CoupleError> {
    Ok(HeatSink {
        name: StableId::new(name)?,
        temperature_k,
    })
}
