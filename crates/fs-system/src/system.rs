//! The system: components, connections, the state store and the step driver.

use core::any::Any;
use core::fmt;
use core::ops::Range;
use std::rc::Rc;

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{
    Carrier, Class, Component, ContractViolation, MassStepAudit, PortLayout, PortRole, PortSpec,
    PortTangent, Proposal, Receipt, STREAM_ENTHALPY_RATE, STREAM_MASS_RATE, audit_mass_step,
    audit_step,
};

use crate::erased::{Node, OwnerError, Plain, Tangent};
use crate::ledger::{BoundaryAccount, ConnectionAccount, Ledger, NodeAccount, SinkAccount, Sum};
use crate::solve::Condensed;
use crate::wiring::{
    AdmissionError, Block, Connection, Driver, HeatSink, NodeId, PortRef, RouteTemperature,
    blocks_in_dependency_order, tolerance_group,
};

/// A prescribed port input: `(t0, dt, values)` fills the value held over the
/// step that starts at `t0` and lasts `dt`.
pub type BoundaryFn = Box<dyn Fn(f64, f64, &mut [f64])>;

/// How many units of first-order rounding a block solve accepts in a
/// residual. The estimate is `eps * sum |d residual / d input| |input|` over
/// the block's own inputs; the factor covers the few operations between an
/// input and a residual and the inputs a block does not solve for.
const ROUNDING_UNITS: f64 = 16.0;

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
        self.need_power_map_roles(a, b)?;
        let (rows, columns) = (self.spec(b).width(), self.spec(a).width());
        if rows.checked_mul(columns) != Some(map.len()) {
            return Err(AdmissionError::MapShape {
                rows,
                columns,
                entries: map.len(),
            });
        }
        if let Some((index, &value)) = map.iter().enumerate().find(|(_, v)| !v.is_finite()) {
            return Err(AdmissionError::MapCoefficient { index, value });
        }
        self.check_free(&[a, b])?;
        let effort_map = self.field_effort_map(a, b, &map)?;
        self.take(&[a, b]);
        self.connections.push(Connection::PowerMap {
            a,
            b,
            map,
            effort_map,
        });
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
                if !taken && spec.role() != PortRole::SignalOut {
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
                Connection::PowerMap {
                    a,
                    b,
                    map,
                    effort_map,
                } => {
                    let columns = self.spec(*a).width();
                    let rows = self.spec(*b).width();
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
                    for (row, index) in input_of(*b).enumerate() {
                        let terms = (0..columns)
                            .filter_map(|column| {
                                let coefficient = map[row * columns + column];
                                (coefficient != 0.0)
                                    .then_some((output_of(*a).start + column, -coefficient))
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
        let block_work = blocks
            .iter()
            .map(|block| {
                // A part without a declared tangent is differenced, and a
                // difference quotient is rounded at the point it is taken
                // at: its slopes are not the same numbers from one step to
                // the next, so factors made from them are not kept.
                let affine = block.nodes.iter().all(|&node| {
                    self.nodes[node].class() == Class::ExactLinear && self.nodes[node].has_tangent()
                });
                BlockWork::plan(block, &drivers, &input_owner, &output_start, affine)
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
                    has_streams: self.nodes[index]
                        .ports()
                        .iter()
                        .any(|port| port.carrier() == Carrier::Stream),
                    ..NodeAccount::default()
                })
                .collect(),
            connections: vec![ConnectionAccount::default(); self.connections.len()],
            boundaries: self
                .boundaries
                .iter()
                .map(|(port, _)| BoundaryAccount {
                    port: self.label(*port),
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
            receipts: vec![Receipt::lossless(0.0, 0.0); count],
            port_energy: vec![0.0; ports],
            port_entropy: vec![0.0; ports],
            port_mass: vec![0.0; ports],
            scratch_outputs: vec![0.0; widest_output],
            scratch_direction: vec![0.0; widest_input],
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
    /// rounding estimate and the configured absolute tolerance. Zero for an
    /// input copied from its driver.
    input_allowance: Vec<f64>,
    receipts: Vec<Receipt>,
    port_energy: Vec<f64>,
    port_entropy: Vec<f64>,
    port_mass: Vec<f64>,
    scratch_outputs: Vec<f64>,
    scratch_direction: Vec<f64>,
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
        let proposal = self.nodes[node].propose(
            self.accepted[node].as_ref(),
            dt_s,
            &self.inputs[inputs],
            self.candidate[node].as_mut(),
            &mut self.outputs[outputs],
            cx,
        );
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
    /// `scratch_outputs`. `Ok(Some(dt))` asks to shorten.
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
        if !matches!(self.nodes[node].class(), Class::ExactLinear | Class::Smooth) {
            return Err(SystemError::Unsupported {
                node: self.names[node].clone(),
                needs: "a port tangent or the non-smooth solve tier",
            });
        }
        // Forward difference of the port map. The candidate slot is free to
        // overwrite: a component without a declared tangent is proposed
        // again, at the accepted inputs, before its block can be accepted.
        let index = inputs.start + local;
        let base = self.inputs[index];
        let scale = base.abs().max(group_scale[self.input_group[index]]);
        let step = f64::EPSILON.sqrt() * if scale > 0.0 { scale } else { 1.0 };
        self.inputs[index] = base + step;
        let proposal = self.nodes[node].propose(
            self.accepted[node].as_ref(),
            dt_s,
            &self.inputs[inputs],
            self.candidate[node].as_mut(),
            &mut self.scratch_outputs[..width],
            cx,
        );
        self.inputs[index] = base;
        match proposal {
            Err(error) => Err(self.refused(node, error, cx)),
            Ok(Proposal::Shorten { dt_s: shorter }) => self.shorter(node, dt_s, shorter),
            Ok(Proposal::Step(_)) => {
                for (column, current) in self.scratch_outputs[..width]
                    .iter_mut()
                    .zip(&self.outputs[outputs])
                {
                    *column = (*column - current) / step;
                }
                Ok(None)
            }
        }
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
    /// rate also the block's largest power). `floor` and `carried` are
    /// scratch of the block's size.
    #[allow(clippy::too_many_arguments)]
    fn carried_beyond(
        &self,
        unknowns: &[usize],
        residual: &[f64],
        slopes: &[(usize, usize, f64)],
        group_scale: &[f64],
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
        for (row, &index) in unknowns.iter().enumerate() {
            let allowed =
                self.config.relative_tolerance * group_scale[self.input_group[index]] + floor[row];
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
        let result = self.solve_lent(&unknowns, &nodes, &mut work, dt_s, cx);
        self.blocks[block].unknowns = unknowns;
        self.blocks[block].nodes = nodes;
        self.block_work[block] = work;
        result
    }

    /// [`Self::solve_block`] on a block's own lists and workspace.
    #[allow(clippy::too_many_lines)]
    fn solve_lent(
        &mut self,
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
        } = work;
        let affine = *affine;
        // Whether this step has already used the factors an affine block
        // kept from an earlier one.
        let mut kept_used = false;
        // Whether this step has assembled the block's Jacobian yet, so that
        // `slopes` are this step's own.
        let mut assembled = false;
        let mut iterations = 0_u32;
        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            for &node in nodes {
                if let Some(shorter) = self.propose(node, dt_s, cx)? {
                    return Ok(Err(shorter));
                }
            }
            // The scale of an unknown is the size of the terms its driver
            // sums, not of their total: a force balance near rest is a small
            // difference of large efforts, and its residual cannot be driven
            // below the roundoff of those efforts.
            group_scale.fill(0.0);
            for (row, &index) in unknowns.iter().enumerate() {
                let driver = &self.drivers[index];
                let target = driver.value(&self.outputs, &self.boundary_values);
                residual[row] = self.inputs[index] - target;
                let group = self.input_group[index];
                group_scale[group] = group_scale[group]
                    .max(self.inputs[index].abs())
                    .max(driver.magnitude(&self.outputs, &self.boundary_values));
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
            for (row, &index) in unknowns.iter().enumerate() {
                // The interface audit is told only this part: it keeps its
                // own relative tolerance, so a loose solver tolerance is
                // still caught there.
                allowance[row] =
                    self.config.absolute_tolerance + ROUNDING_UNITS * f64::EPSILON * rounding[row];
                if self.input_is_energy_rate[index] {
                    allowance[row] += ROUNDING_UNITS * f64::EPSILON * power_scale;
                }
                let allowed = self.config.relative_tolerance * group_scale[self.input_group[index]]
                    + allowance[row];
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
                        group_scale,
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
                    if let Some(shorter) =
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
                    group_scale,
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
                        self.input_allowance[index] = told;
                    }
                    return Ok(Ok(iterations));
                };
                if iterations >= self.config.max_iterations {
                    return Err(self.not_converged(nodes, unknowns, residual, iterations, open));
                }
            }
            iterations += 1;
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
                    &self.inputs[inputs],
                    &self.outputs[outputs],
                    old_mass,
                    new_mass,
                )
                .map_err(|violation| SystemError::Contract {
                    node: self.names[node].clone(),
                    violation,
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
        }
        for (account, &(energy, entropy, mass)) in
            self.ledger.connections.iter_mut().zip(&connection_rows)
        {
            account.energy_defect_j.add(energy);
            account.entropy_generated_j_per_k.add(entropy);
            account.mass_defect_kg.add(mass);
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
        }
        Ok(StepReport {
            dt_s,
            time_s: self.time.value(),
            iterations,
            shortenings,
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
