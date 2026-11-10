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
    /// Two connected ports carry different things or have different widths.
    Mismatch {
        /// First port, as `node.port`.
        a: String,
        /// Second port, as `node.port`.
        b: String,
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
                write!(f, "ports `{a}` and `{b}` differ in carrier or width")
            }
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
        }
    }
}

impl std::error::Error for AdmissionError {}

/// A lossless or routing relation between ports.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Connection {
    /// The setter returns the effort every member receives; the setter
    /// receives minus the sum of the members' flows.
    CommonEffort {
        setter: PortRef,
        members: Vec<PortRef>,
    },
    /// The setter returns the flow; every member receives minus that flow and
    /// the setter receives the sum of the members' efforts.
    CommonFlow {
        setter: PortRef,
        members: Vec<PortRef>,
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
                out.extend_from_slice(members);
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
            Self::Linear(terms) => {
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
            Self::Linear(terms) => {
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

    /// Derivative of the value with respect to each output it reads.
    pub(crate) fn gradient(&self, outputs: &[f64], out: &mut Vec<(usize, f64)>) {
        out.clear();
        match self {
            Self::Linear(terms) => out.extend_from_slice(terms),
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
            Self::Linear(terms) | Self::EquivalentTemperature { terms, .. } => {
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
/// the largest of their kind.
pub(crate) fn tolerance_group(groups: &mut Vec<(Carrier, PortRole)>, port: &PortSpec) -> usize {
    let key = (port.carrier(), port.role());
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
