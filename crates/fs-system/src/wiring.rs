//! Connections between ports, their admission, and the evaluation plan.
//!
//! A connection is data: the kernel generates its equations. Every scalar
//! input of every component gets exactly one [`Driver`], which says how the
//! input is computed from outputs, from a boundary, or from a constant. The
//! drivers define which component depends on which, and the strongly
//! connected blocks of that dependency graph, in order, are the plan.

use core::fmt;

use fs_ports::StableId;
use fs_ports::component::{Carrier, ContractViolation, PortRole, PortSpec};

/// A component inside a system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(pub(crate) usize);

impl NodeId {
    /// Position of the component in the system and in the ledger's accounts.
    #[must_use]
    pub fn index(self) -> usize {
        self.0
    }
}

/// One port of one component inside a system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PortRef {
    pub(crate) node: usize,
    pub(crate) port: usize,
}

/// Which temperature the source of a heat route receives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RouteTemperature {
    /// The entropy-equivalent temperature of the receivers,
    /// `sum(f_i) / sum(f_i / T_i)`. The route then generates no entropy, and
    /// the source depends on the receivers.
    SinkEquivalent,
    /// A fixed temperature [K], for a source whose law does not use it. The
    /// route is then one-way: the receivers depend on the source only.
    Fixed(f64),
}

/// A named place outside the system that receives heat nobody else was given.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatSink {
    /// Name under which the ledger reports the heat.
    pub name: StableId,
    /// Temperature [K] at which the heat is deposited.
    pub temperature_k: f64,
}

/// Why a system description was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum AdmissionError {
    /// Two components share a name.
    DuplicateNode {
        /// The repeated name.
        name: String,
    },
    /// A component has no port with this identifier.
    UnknownPort {
        /// Component name.
        node: String,
        /// Port identifier asked for.
        port: String,
    },
    /// A port list is not a valid layout.
    Layout {
        /// Component name.
        node: String,
        /// What is wrong with its ports.
        violation: ContractViolation,
    },
    /// A port has the wrong role for its place in a connection.
    Role {
        /// The port, as `node.port`.
        port: String,
        /// What the connection needs there.
        needs: &'static str,
    },
    /// Connected ports differ in carrier, width or power-pairing metric.
    Mismatch {
        /// First port, as `node.port`.
        a: String,
        /// Second port, as `node.port`.
        b: String,
    },
    /// Two heat routes name the same sink at different temperatures.
    SinkTemperature {
        /// The sink's name.
        name: String,
        /// Temperature [K] of the first declaration.
        declared_k: f64,
        /// Temperature [K] of the refused one.
        repeated_k: f64,
    },
    /// A port is used by more than one connection or boundary.
    AlreadyConnected {
        /// The port, as `node.port`.
        port: String,
    },
    /// A port that exchanges power or needs an input is connected to nothing.
    Unconnected {
        /// The port, as `node.port`.
        port: String,
    },
    /// Heat-route fractions are not in `(0, 1]` or sum to more than one.
    Fractions {
        /// Sum of the declared fractions.
        sum: f64,
    },
    /// Heat-route fractions sum to less than one and no sink was declared.
    SinkRequired {
        /// Fraction that would go nowhere.
        unallocated: f64,
    },
    /// A temperature is not positive and finite.
    Temperature {
        /// The rejected value [K].
        kelvin: f64,
    },
    /// A connection has no member.
    Empty,
    /// A common-flow junction of heat ports with more than one member. Heat
    /// rate is power itself, so giving the same heat rate to several members
    /// would multiply the energy.
    HeatCommonFlow {
        /// The setter port, as `node.port`.
        port: String,
    },
    /// A common-flow junction of stream ports with more than one member. A
    /// stream's rates are what crosses, so several storages cannot each
    /// receive the whole of one flow. Streams meet at a storage.
    StreamCommonFlow {
        /// The setter port, as `node.port`.
        port: String,
    },
    /// An ideal transformer's ratio is zero or not finite.
    TransformerRatio {
        /// The rejected signed ratio.
        ratio: f64,
    },
    /// A dense power map does not have one entry per row and column.
    MapShape {
        /// Width of the flow-receiving port.
        rows: usize,
        /// Width of the effort-receiving port.
        columns: usize,
        /// Number of coefficients supplied.
        entries: usize,
    },
    /// A power-map coefficient is not finite.
    MapCoefficient {
        /// Row-major index of the rejected coefficient.
        index: usize,
        /// The rejected value.
        value: f64,
    },
    /// A field map's weighted effort adjoint could not be represented.
    MapAdjoint {
        /// Effort-receiving port, as `node.port`.
        a: String,
        /// Flow-receiving port, as `node.port`.
        b: String,
        /// The rejected metric operation or non-finite intermediate.
        violation: ContractViolation,
    },
    /// Two components on a load line are in one loop. The solve finds one
    /// line per loop.
    TwoLines {
        /// The first of them, by name.
        first: String,
        /// The second.
        second: String,
    },
    /// A component's load line is not over power ports that receive a flow,
    /// listed in increasing order.
    LinePorts {
        /// Component name.
        node: String,
    },
}

impl fmt::Display for AdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateNode { name } => write!(f, "component name `{name}` is used twice"),
            Self::UnknownPort { node, port } => {
                write!(f, "component `{node}` has no port `{port}`")
            }
            Self::Layout { node, violation } => write!(f, "component `{node}`: {violation}"),
            Self::Role { port, needs } => write!(f, "port `{port}` must be {needs} here"),
            Self::Mismatch { a, b } => {
                write!(
                    f,
                    "ports `{a}` and `{b}` differ in carrier, width or power metric"
                )
            }
            Self::SinkTemperature {
                name,
                declared_k,
                repeated_k,
            } => write!(
                f,
                "sink `{name}` was declared at {declared_k} K and is named again at {repeated_k} K"
            ),
            Self::AlreadyConnected { port } => write!(f, "port `{port}` is already connected"),
            Self::Unconnected { port } => write!(f, "port `{port}` is connected to nothing"),
            Self::Fractions { sum } => {
                write!(
                    f,
                    "heat-route fractions must lie in (0, 1] and sum to at most 1 (sum {sum})"
                )
            }
            Self::SinkRequired { unallocated } => write!(
                f,
                "heat-route fractions leave {unallocated} unallocated and no sink is declared"
            ),
            Self::Temperature { kelvin } => {
                write!(f, "temperature {kelvin} K is not positive and finite")
            }
            Self::Empty => write!(f, "a connection needs at least one member"),
            Self::HeatCommonFlow { port } => write!(
                f,
                "heat port `{port}` cannot set a common flow for several members: heat ports \
                 share a temperature (common effort) or use a heat route"
            ),
            Self::StreamCommonFlow { port } => write!(
                f,
                "stream port `{port}` cannot set a common flow for several members: streams \
                 meet at a storage, which sets what they share (common effort)"
            ),
            Self::TransformerRatio { ratio } => {
                write!(
                    f,
                    "transformer ratio must be finite and nonzero (got {ratio})"
                )
            }
            Self::MapShape {
                rows,
                columns,
                entries,
            } => write!(
                f,
                "power map needs a {rows} by {columns} row-major matrix (got {entries} entries)"
            ),
            Self::MapCoefficient { index, value } => {
                write!(
                    f,
                    "power-map coefficient {index} must be finite (got {value})"
                )
            }
            Self::MapAdjoint { a, b, violation } => {
                write!(
                    f,
                    "power map `{a}` to `{b}` has no representable effort adjoint: {violation}"
                )
            }
            Self::TwoLines { first, second } => write!(
                f,
                "components `{first}` and `{second}` are both on a load line in one loop; the \
                 solve finds one line per loop"
            ),
            Self::LinePorts { node } => write!(
                f,
                "the load line of component `{node}` must span power ports that receive a flow, \
                 in increasing order"
            ),
        }
    }
}

impl std::error::Error for AdmissionError {}

/// A lossless or routing relation between ports.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Connection {
    /// The setter returns the effort; member `i` receives `ratio_i` times it
    /// and the setter receives minus the sum of `ratio_i` times the members'
    /// flows. Every ratio of a plain junction is one.
    CommonEffort {
        setter: PortRef,
        members: Vec<(PortRef, f64)>,
    },
    /// The setter returns the flow; member `i` receives minus `ratio_i` times
    /// it and the setter receives the sum of `ratio_i` times the members'
    /// efforts. Every ratio of a plain junction is one.
    CommonFlow {
        setter: PortRef,
        members: Vec<(PortRef, f64)>,
    },
    /// `b` receives `-ratio * flow_a`; `a` receives `ratio * effort_b`.
    /// A signed scalar ratio acts independently on each coefficient.
    Transformer { a: PortRef, b: PortRef, ratio: f64 },
    /// `b` receives `-P flow_a`; `a` receives the metric adjoint times
    /// `effort_b`. `map` has `b.width()` rows and `a.width()` columns.
    PowerMap {
        a: PortRef,
        b: PortRef,
        map: Vec<f64>,
        /// Row-major `M_a^-1 P^T M_b`, with `a.width()` rows. `None`
        /// retains the original shared coefficients for ordinary ports.
        effort_map: Option<Vec<f64>>,
    },
    /// Heat leaving `source` is delivered to `sinks` by fraction; the rest
    /// goes to the declared sink outside the system.
    HeatRoute {
        source: PortRef,
        sinks: Vec<(PortRef, f64)>,
        /// `(index into the system's sink list, fraction, temperature)`.
        remainder: Option<(usize, f64, f64)>,
        temperature: RouteTemperature,
    },
    /// A one-way copy of a signal.
    Signal { from: PortRef, to: PortRef },
}

impl Connection {
    /// Every port that exchanges energy through this connection.
    pub(crate) fn energy_members(&self, out: &mut Vec<PortRef>) {
        out.clear();
        match self {
            Self::CommonEffort { setter, members } | Self::CommonFlow { setter, members } => {
                out.push(*setter);
                out.extend(members.iter().map(|(port, _)| *port));
            }
            Self::Transformer { a, b, .. } | Self::PowerMap { a, b, .. } => {
                out.push(*a);
                out.push(*b);
            }
            Self::HeatRoute { source, sinks, .. } => {
                out.push(*source);
                out.extend(sinks.iter().map(|(port, _)| *port));
            }
            Self::Signal { .. } => {}
        }
    }
}

/// How one scalar input is obtained.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Driver {
    /// Nothing drives it yet (refused at sealing).
    Unset,
    /// `sum(coefficient * outputs[index])`.
    Linear(Vec<(usize, f64)>),
    /// The same sum, for a flow on a load line. The last `own` terms read
    /// the line component's own efforts, and their coefficients are the
    /// line's compliance, which the solve sets anew in every step; with
    /// them the value is the released flow, the flow the port would have
    /// if the component returned no effort.
    OnLine {
        terms: Vec<(usize, f64)>,
        own: usize,
    },
    /// Entropy-equivalent temperature of a heat route's receivers:
    /// `total / (sum(fraction / outputs[index]) + fixed)`.
    EquivalentTemperature {
        terms: Vec<(usize, f64)>,
        /// `fraction / temperature` of the declared sink, or zero.
        fixed: f64,
        total: f64,
    },
    /// Component `component` of boundary `boundary`.
    Boundary { boundary: usize, component: usize },
    /// A constant.
    Fixed(f64),
}

impl Driver {
    /// Value of the input given all current outputs and boundary values.
    pub(crate) fn value(&self, outputs: &[f64], boundaries: &[Vec<f64>]) -> f64 {
        match self {
            Self::Unset => f64::NAN,
            Self::Linear(terms) | Self::OnLine { terms, .. } => {
                let mut sum = 0.0;
                for &(index, coefficient) in terms {
                    sum += coefficient * outputs[index];
                }
                sum
            }
            Self::EquivalentTemperature {
                terms,
                fixed,
                total,
            } => {
                let mut inverse = *fixed;
                for &(index, fraction) in terms {
                    inverse += fraction / outputs[index];
                }
                total / inverse
            }
            Self::Boundary {
                boundary,
                component,
            } => boundaries[*boundary][*component],
            Self::Fixed(value) => *value,
        }
    }

    /// Size of the terms that make up the value: the sum of the magnitudes a
    /// linear driver adds, or the magnitude of the value itself. A convergence
    /// test scaled by the value alone cannot be met when the terms cancel.
    pub(crate) fn magnitude(&self, outputs: &[f64], boundaries: &[Vec<f64>]) -> f64 {
        match self {
            Self::Linear(terms) | Self::OnLine { terms, .. } => {
                let mut sum = 0.0;
                for &(index, coefficient) in terms {
                    sum += (coefficient * outputs[index]).abs();
                }
                sum
            }
            Self::Unset
            | Self::EquivalentTemperature { .. }
            | Self::Boundary { .. }
            | Self::Fixed(_) => self.value(outputs, boundaries).abs(),
        }
    }

    /// [`Self::magnitude`] without a load line's own terms: the size of the
    /// flow the other components give a port on a line, which is what
    /// crosses there together with what the line component takes away.
    pub(crate) fn magnitude_of_others(&self, outputs: &[f64], boundaries: &[Vec<f64>]) -> f64 {
        match self {
            Self::OnLine { terms, own } => {
                let mut sum = 0.0;
                for &(index, coefficient) in &terms[..terms.len() - own] {
                    sum += (coefficient * outputs[index]).abs();
                }
                sum
            }
            _ => self.magnitude(outputs, boundaries),
        }
    }

    /// Derivative of the value with respect to each output it reads.
    pub(crate) fn gradient(&self, outputs: &[f64], out: &mut Vec<(usize, f64)>) {
        out.clear();
        match self {
            Self::Linear(terms) | Self::OnLine { terms, .. } => out.extend_from_slice(terms),
            Self::EquivalentTemperature {
                terms,
                fixed,
                total,
            } => {
                let mut inverse = *fixed;
                for &(index, fraction) in terms {
                    inverse += fraction / outputs[index];
                }
                for &(index, fraction) in terms {
                    let t = outputs[index];
                    out.push((index, total * fraction / (t * t * inverse * inverse)));
                }
            }
            Self::Unset | Self::Boundary { .. } | Self::Fixed(_) => {}
        }
    }

    /// Indices of the outputs this driver reads.
    pub(crate) fn sources(&self, out: &mut Vec<usize>) {
        match self {
            Self::Linear(terms)
            | Self::OnLine { terms, .. }
            | Self::EquivalentTemperature { terms, .. } => {
                out.extend(terms.iter().map(|(index, _)| *index));
            }
            Self::Unset | Self::Boundary { .. } | Self::Fixed(_) => {}
        }
    }
}

/// One step of the plan: components solved together.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Block {
    /// Components, in index order.
    pub(crate) nodes: Vec<usize>,
    /// Global input indices that depend on outputs of this block.
    pub(crate) unknowns: Vec<usize>,
    /// The other inputs of its components, which earlier blocks, boundaries
    /// and constants supply, in component and then port order.
    pub(crate) known: Vec<usize>,
}

impl Block {
    /// Whether the block needs an iterative solve.
    pub(crate) fn looped(&self) -> bool {
        !self.unknowns.is_empty()
    }
}

/// Strongly connected components of `successors`, each sorted, in an order
/// where every block comes after the blocks it depends on.
///
/// `successors[i]` lists the nodes that depend on node `i`. Iterative Tarjan,
/// visiting nodes and edges in index order, so the result is deterministic.
pub(crate) fn blocks_in_dependency_order(successors: &[Vec<usize>]) -> Vec<Vec<usize>> {
    const UNSEEN: usize = usize::MAX;
    let n = successors.len();
    let mut index = vec![UNSEEN; n];
    let mut low = vec![0_usize; n];
    let mut on_stack = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut work: Vec<(usize, usize)> = Vec::new();
    let mut next = 0_usize;
    let mut found: Vec<Vec<usize>> = Vec::new();
    for root in 0..n {
        if index[root] != UNSEEN {
            continue;
        }
        work.push((root, 0));
        while let Some(&(node, edge)) = work.last() {
            if edge == 0 {
                index[node] = next;
                low[node] = next;
                next += 1;
                stack.push(node);
                on_stack[node] = true;
            }
            if let Some(&successor) = successors[node].get(edge) {
                if let Some(top) = work.last_mut() {
                    top.1 += 1;
                }
                if index[successor] == UNSEEN {
                    work.push((successor, 0));
                } else if on_stack[successor] {
                    low[node] = low[node].min(index[successor]);
                }
                continue;
            }
            work.pop();
            if let Some(&(parent, _)) = work.last() {
                low[parent] = low[parent].min(low[node]);
            }
            if low[node] == index[node] {
                let mut block = Vec::new();
                while let Some(member) = stack.pop() {
                    on_stack[member] = false;
                    block.push(member);
                    if member == node {
                        break;
                    }
                }
                block.sort_unstable();
                found.push(block);
            }
        }
    }
    // Tarjan emits a block only after every block reachable from it, that is,
    // after its dependents. Reversing puts dependencies first.
    found.reverse();
    found
}

/// A tolerance group: unknowns of the same carrier and side are compared with
/// the largest of their kind. `offset` is the coefficient's position in its
/// port; it tells the quantities of a stream port apart and is ignored for
/// the other carriers, whose coefficients are all one quantity.
pub(crate) fn tolerance_group(
    groups: &mut Vec<(Carrier, PortRole, usize)>,
    port: &PortSpec,
    offset: usize,
) -> usize {
    let position = if port.carrier() == Carrier::Stream {
        offset
    } else {
        0
    };
    let key = (port.carrier(), port.role(), position);
    if let Some(position) = groups.iter().position(|known| *known == key) {
        position
    } else {
        groups.push(key);
        groups.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use super::{Driver, blocks_in_dependency_order};

    #[test]
    fn blocks_follow_dependencies_and_group_cycles() {
        // 0 -> 1 <-> 2 -> 3, and 4 alone.
        let successors = vec![vec![1], vec![2], vec![1, 3], vec![], vec![]];
        let blocks = blocks_in_dependency_order(&successors);
        let position = |node: usize| blocks.iter().position(|b| b.contains(&node)).unwrap();
        assert!(blocks.contains(&vec![1, 2]));
        assert!(position(0) < position(1));
        assert!(position(2) < position(3));
        assert_eq!(blocks.iter().map(Vec::len).sum::<usize>(), 5);
        // Same input, same plan.
        assert_eq!(blocks, blocks_in_dependency_order(&successors));
    }

    #[test]
    fn a_long_chain_does_not_recurse() {
        let n = 200_000;
        let successors: Vec<Vec<usize>> = (0..n)
            .map(|i| if i + 1 < n { vec![i + 1] } else { vec![] })
            .collect();
        let blocks = blocks_in_dependency_order(&successors);
        assert_eq!(blocks.len(), n);
        assert_eq!(blocks[0], vec![0]);
        assert_eq!(blocks[n - 1], vec![n - 1]);
    }

    #[test]
    fn equivalent_temperature_gradient_matches_a_difference() {
        let driver = Driver::EquivalentTemperature {
            terms: vec![(0, 0.6), (1, 0.3)],
            fixed: 0.1 / 293.0,
            total: 1.0,
        };
        let outputs = [310.0, 345.0];
        let mut gradient = Vec::new();
        driver.gradient(&outputs, &mut gradient);
        for &(index, slope) in &gradient {
            let mut up = outputs;
            let mut down = outputs;
            up[index] += 1.0e-3;
            down[index] -= 1.0e-3;
            let difference = (driver.value(&up, &[]) - driver.value(&down, &[])) / 2.0e-3;
            assert!((difference - slope).abs() <= 1.0e-8 * slope.abs());
        }
        // Between the coldest and the hottest receiver.
        let t = driver.value(&outputs, &[]);
        assert!(t > 293.0 && t < 345.0);
    }
}
