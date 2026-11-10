//! The system: components, connections, the state store and the step driver.

use core::any::Any;
use core::fmt;
use core::ops::Range;

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{
    Carrier, Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, PortTangent,
    Proposal, Receipt, audit_step,
};

use crate::erased::{Node, OwnerError, Plain, Tangent};
use crate::ledger::{BoundaryAccount, ConnectionAccount, Ledger, NodeAccount, SinkAccount, Sum};
use crate::solve::solve_dense;
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

/// Limits of the step driver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveConfig {
    /// Interface residual allowed, relative to the larger of the two sides
    /// and to the largest unknown of the same kind in the block.
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
    /// `stored_energy_j`, during the audit.
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
        /// The input with the largest scaled residual, as `node.port`.
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
            Ok(())
        } else {
            Err(AdmissionError::Mismatch {
                a: self.label(a),
                b: self.label(b),
            })
        }
    }

    fn junction(
        &mut self,
        setter: PortRef,
        members: &[PortRef],
        setter_role: PortRole,
        member_role: PortRole,
    ) -> Result<(), AdmissionError> {
        if members.is_empty() {
            return Err(AdmissionError::Empty);
        }
        let (setter_needs, member_needs) = if setter_role == PortRole::FlowIn {
            ("a port that receives flow", "a port that receives effort")
        } else {
            ("a port that receives effort", "a port that receives flow")
        };
        self.need_role(setter, setter_role, setter_needs)?;
        for &member in members {
            self.need_role(member, member_role, member_needs)?;
            self.need_match(setter, member)?;
        }
        // Heat rate is power itself, so handing the same heat rate to several
        // members would multiply the energy. Heat ports share a temperature
        // (common effort) or go through a heat route.
        if setter_role == PortRole::EffortIn
            && members.len() > 1
            && self.spec(setter).carrier() == Carrier::Heat
        {
            return Err(AdmissionError::HeatCommonFlow {
                port: self.label(setter),
            });
        }
        let mut all = Vec::with_capacity(members.len() + 1);
        all.push(setter);
        all.extend_from_slice(members);
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
        self.junction(setter, members, PortRole::FlowIn, PortRole::EffortIn)?;
        self.connections.push(Connection::CommonEffort {
            setter,
            members: members.to_vec(),
        });
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
        self.junction(setter, members, PortRole::EffortIn, PortRole::FlowIn)?;
        self.connections.push(Connection::CommonFlow {
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
                    // the setter, the setter gets the members' sum.
                    let (to_member, to_setter) = if common_effort {
                        (1.0, -1.0)
                    } else {
                        (-1.0, 1.0)
                    };
                    let setter_out = output_of(*setter);
                    for (offset, index) in input_of(*setter).enumerate() {
                        let terms = members
                            .iter()
                            .map(|member| (output_of(*member).start + offset, to_setter))
                            .collect();
                        drivers[index] = Driver::Linear(terms);
                    }
                    for member in members {
                        for (offset, index) in input_of(*member).enumerate() {
                            drivers[index] =
                                Driver::Linear(vec![(setter_out.start + offset, to_member)]);
                        }
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

        // Tolerance group and label of every input. Power and heat inputs of
        // one carrier and side are compared with the largest of their kind in
        // the block. A signal has no kind, so each scalar signal input is a
        // group of its own and is judged against its own size only.
        let mut groups = Vec::new();
        let mut signal_inputs = Vec::new();
        let mut input_group = vec![0_usize; inputs];
        let mut input_label = vec![String::new(); inputs];
        let mut input_is_temperature = vec![false; inputs];
        for node in 0..count {
            for (port, spec) in self.nodes[node].ports().iter().enumerate() {
                let reference = PortRef { node, port };
                for index in input_of(reference) {
                    if spec.carrier() == Carrier::Signal {
                        signal_inputs.push(index);
                    } else {
                        input_group[index] = tolerance_group(&mut groups, spec);
                    }
                    input_label[index] = self.label(reference);
                    input_is_temperature[index] =
                        spec.carrier() == Carrier::Heat && spec.role() == PortRole::EffortIn;
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
                .map(|name| NodeAccount {
                    name: name.clone(),
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
                let Connection::HeatRoute { sinks: theirs, .. } = &self.connections[earlier]
                else {
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
            input_start,
            output_start,
            port_start,
            input_owner,
            input_group,
            input_label,
            input_is_temperature,
            group_count,
            inputs: vec![0.0; inputs],
            outputs: vec![0.0; outputs],
            accepted_inputs: vec![0.0; inputs],
            accepted_outputs: vec![0.0; outputs],
            input_allowance: vec![0.0; inputs],
            receipts: vec![Receipt::lossless(0.0, 0.0); count],
            port_energy: vec![0.0; ports],
            port_entropy: vec![0.0; ports],
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

/// Components joined by connections, stepped as one.
///
/// The system owns every state in two generations, accepted and candidate. A
/// step fills the candidates, audits and admits them, and only then publishes
/// all of them by swapping; a failed step leaves every accepted state, the
/// clock and the ledger untouched.
pub struct System {
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
    input_start: Vec<usize>,
    output_start: Vec<usize>,
    port_start: Vec<usize>,
    input_owner: Vec<usize>,
    input_group: Vec<usize>,
    input_label: Vec<String>,
    input_is_temperature: Vec<bool>,
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
    scratch_outputs: Vec<f64>,
    scratch_direction: Vec<f64>,
    time: Sum,
    steps: u64,
    ledger: Ledger,
    config: SolveConfig,
}

impl System {
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

    fn block_names(&self, block: usize) -> Vec<String> {
        self.blocks[block]
            .nodes
            .iter()
            .map(|&node| self.names[node].clone())
            .collect()
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
        // overwrite: the block is evaluated again before it can be accepted.
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

    /// Solve one looped block by Newton iteration on its interface inputs.
    #[allow(clippy::too_many_lines)]
    fn solve_block(
        &mut self,
        block: usize,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        let unknowns = self.blocks[block].unknowns.clone();
        let nodes = self.blocks[block].nodes.clone();
        let m = unknowns.len();
        for &index in &unknowns {
            let start = self.accepted_inputs[index];
            if self.input_is_temperature[index] && !(start.is_finite() && start > 0.0) {
                return Err(SystemError::GuessRequired {
                    port: self.input_label[index].clone(),
                });
            }
            self.inputs[index] = start;
        }
        let mut residual = vec![0.0; m];
        let mut jacobian = vec![0.0; m * m];
        // First-order size of the rounding in each residual: how strongly its
        // driver's sources depend on the block's inputs, times the size of
        // those inputs. Known after the first tangent assembly.
        let mut rounding = vec![0.0_f64; m];
        let mut allowance = vec![0.0_f64; m];
        let mut group_scale = vec![0.0_f64; self.group_count];
        let mut gradient = Vec::new();
        let mut iterations = 0_u32;
        loop {
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            for &node in &nodes {
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
            let mut worst = (0_usize, 0.0_f64);
            let mut converged = true;
            for (row, &index) in unknowns.iter().enumerate() {
                // The interface audit is told only this part: it keeps its
                // own relative tolerance, so a loose solver tolerance is
                // still caught there.
                allowance[row] = self.config.absolute_tolerance
                    + ROUNDING_UNITS * f64::EPSILON * rounding[row];
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
            if converged {
                for (row, &index) in unknowns.iter().enumerate() {
                    self.input_allowance[index] = allowance[row];
                }
                return Ok(Ok(iterations));
            }
            if iterations >= self.config.max_iterations {
                return Err(SystemError::NotConverged {
                    nodes: self.block_names(block),
                    iterations,
                    worst: self.input_label[unknowns[worst.0]].clone(),
                    residual: residual[worst.0],
                });
            }
            iterations += 1;

            // J = I - (d driver / d outputs) (d outputs / d inputs).
            rounding.fill(0.0);
            jacobian.fill(0.0);
            for diagonal in 0..m {
                jacobian[diagonal * m + diagonal] = 1.0;
            }
            for (column, &index) in unknowns.iter().enumerate() {
                let node = self.input_owner[index];
                let local = index - self.input_start[node];
                if let Some(shorter) = self.tangent_column(node, local, dt_s, &group_scale, cx)? {
                    return Ok(Err(shorter));
                }
                let outputs = self.node_outputs(node);
                for (row, &target) in unknowns.iter().enumerate() {
                    self.drivers[target].gradient(&self.outputs, &mut gradient);
                    for &(source, slope) in &gradient {
                        if outputs.contains(&source) {
                            let through = slope * self.scratch_outputs[source - outputs.start];
                            jacobian[row * m + column] -= through;
                            rounding[row] += (through * self.inputs[index]).abs();
                        }
                    }
                }
            }
            for value in &mut residual {
                *value = -*value;
            }
            if !solve_dense(&mut jacobian, &mut residual, m) {
                return Err(SystemError::Singular {
                    nodes: self.block_names(block),
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
            for position in 0..self.blocks[block].nodes.len() {
                let node = self.blocks[block].nodes[position];
                for index in self.node_inputs(node) {
                    if !self.blocks[block].unknowns.contains(&index) {
                        self.inputs[index] =
                            self.drivers[index].value(&self.outputs, &self.boundary_values);
                    }
                }
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
    fn audit(&mut self, dt_s: f64, defects: &mut [f64]) -> Result<(), SystemError> {
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

        let mut defects = vec![0.0; self.nodes.len()];
        self.audit(dt_s, &mut defects)?;

        // Connection and boundary rows of this step.
        let mut members = Vec::new();
        let mut connection_rows = Vec::with_capacity(self.connections.len());
        let mut sink_heat = vec![0.0; self.ledger.sinks.len()];
        let (mut generated, mut entropy_scale) = (0.0_f64, 0.0_f64);
        for receipt in &self.receipts {
            generated += receipt.entropy_generated_j_per_k;
            entropy_scale += receipt.entropy_generated_j_per_k;
        }
        let mut group_members = Vec::new();
        for (index, group) in self.audit_groups.iter().enumerate() {
            // A heat route that shares a receiver with an earlier one is
            // audited in that one's row.
            if group.is_empty() {
                connection_rows.push((0.0, 0.0));
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
            let (mut energy, mut magnitude, mut entropy) = (0.0_f64, 0.0_f64, 0.0_f64);
            for port in &members {
                let global = self.port_start[port.node] + port.port;
                energy += self.port_energy[global];
                magnitude += self.port_energy[global].abs();
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
            // heat rate itself.
            let mut solve_slack_j = 0.0_f64;
            for port in &members {
                let spec = &self.nodes[port.node].ports()[port.port];
                let weight_one = match (spec.carrier(), spec.role()) {
                    (Carrier::Heat, PortRole::FlowIn) => Some(1.0),
                    (Carrier::Power(_), _) => None,
                    _ => Some(0.0),
                };
                let outputs = self.output_range(*port);
                for (offset, index) in self.input_range(*port).enumerate() {
                    let weight = weight_one
                        .unwrap_or_else(|| self.outputs[outputs.start + offset].abs());
                    solve_slack_j += dt_s * self.input_allowance[index] * weight;
                }
            }
            let tolerance_j = self.config.interface_energy_tolerance * magnitude
                + solve_slack_j
                + f64::MIN_POSITIVE;
            if !(energy.abs() <= tolerance_j) {
                return Err(SystemError::Interface {
                    connection: index,
                    defect_j: energy,
                    tolerance_j,
                });
            }
            generated += entropy;
            connection_rows.push((energy, entropy));
        }
        if !(generated >= -1.0e-9 * entropy_scale) {
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
        }
        for (account, (energy, entropy)) in self.ledger.connections.iter_mut().zip(connection_rows)
        {
            account.energy_defect_j.add(energy);
            account.entropy_generated_j_per_k.add(entropy);
        }
        for (account, heat) in self.ledger.sinks.iter_mut().zip(sink_heat) {
            account.heat_j.add(heat);
        }
        for (account, (port, _)) in self.ledger.boundaries.iter_mut().zip(&self.boundaries) {
            let global = self.port_start[port.node] + port.port;
            account.energy_in_j.add(self.port_energy[global]);
            account.entropy_in_j_per_k.add(self.port_entropy[global]);
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
