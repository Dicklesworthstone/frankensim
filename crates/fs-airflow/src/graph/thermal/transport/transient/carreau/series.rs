//! Explicit zero-volume continuations between original finite pipe sections.
//!
//! Only thermal topology is condensed. Original hydraulic nodes, branches,
//! local rheology, radial histories and work accounts retain their identities.

use super::{
    CarreauGraphError, CarreauPipeGraph, Cx, Density, FixedPressure, FlowInjection,
    ImplicitTransportError, thermal_checkpoint,
};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

#[derive(Clone)]
pub(super) struct SeriesInventory {
    pub(super) nodes: Vec<Option<CellId>>,
    pub(super) junctions: Vec<CellId>,
    pub(super) connections: Vec<Connection>,
    splices: Vec<Splice>,
}

#[derive(Clone)]
pub(super) struct Connection {
    pub(super) from: usize,
    pub(super) to: usize,
    pub(super) branch: usize,
    pub(super) sign: f64,
}

#[derive(Clone, Copy)]
struct Incident {
    branch: usize,
    /// Positive signed flow directed into this hydraulic node.
    inward_sign: f64,
}

#[derive(Clone)]
struct Splice {
    incidents: [Incident; 2],
}

impl SeriesInventory {
    /// Parent has admitted the full branch roster, real pipe geometry and all
    /// coincident hydraulic endpoints. No fictional cell is used for admission.
    pub(super) fn new(
        cx: &Cx<'_>,
        graph: &CarreauPipeGraph,
        nodes: &[Option<CellId>],
        pipe_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Self, ImplicitTransportError<CarreauGraphError>> {
        let junctions: Vec<_> = nodes.iter().copied().flatten().collect();
        let mut slots = Vec::with_capacity(nodes.len());
        let mut slot = 0;
        for node in nodes {
            thermal_checkpoint(cx)?;
            slots.push(node.map(|_| {
                let result = slot;
                slot += 1;
                result
            }));
        }
        let mut incidents = vec![Vec::new(); nodes.len()];
        for (index, branch) in graph.branches().iter().enumerate() {
            thermal_checkpoint(cx)?;
            for (node, sign) in [(branch.from, -1.0), (branch.to, 1.0)] {
                if nodes[node].is_none() {
                    incidents[node].push(Incident {
                        branch: index,
                        inward_sign: sign,
                    });
                    if incidents[node].len() > 2 {
                        return Err(ImplicitTransportError::Invalid(
                            "zero-volume Carreau splice requires exactly two incident pipes",
                        ));
                    }
                }
            }
        }
        let mut splices = Vec::new();
        for (node, binding) in nodes.iter().enumerate() {
            thermal_checkpoint(cx)?;
            if binding.is_some() {
                continue;
            }
            let [first, second] = incidents[node].as_slice() else {
                return Err(ImplicitTransportError::Invalid(
                    "zero-volume Carreau splice requires exactly two incident pipes",
                ));
            };
            let a = &graph.branches()[first.branch];
            let b = &graph.branches()[second.branch];
            if first.branch == second.branch
                || a.pipe.radius_m() != b.pipe.radius_m()
                || a.pipe.law() != b.pipe.law()
            {
                return Err(ImplicitTransportError::Invalid(
                    "Carreau series splice requires distinct same-radius original-fluid sections",
                ));
            }
            let direction = |incident: &Incident| {
                let geometry = thermal.cells()[&pipe_cells[incident.branch]].geometry;
                let delta = std::array::from_fn::<_, 3, _>(|i| {
                    -incident.inward_sign * (geometry.end_m[i] - geometry.start_m[i])
                });
                let length = fs_math::det::hypot(fs_math::det::hypot(delta[0], delta[1]), delta[2]);
                delta.map(|value| value / length)
            };
            let first_direction = direction(first);
            let second_direction = direction(second);
            if (0..3).any(|axis| {
                !(first_direction[axis] + second_direction[axis]).is_finite()
                    || (first_direction[axis] + second_direction[axis]).abs() > 64.0 * f64::EPSILON
            }) {
                return Err(ImplicitTransportError::Invalid(
                    "Carreau series splice requires straight contiguous native pipe centrelines",
                ));
            }
            splices.push(Splice {
                incidents: [*first, *second],
            });
        }
        let mut connections = Vec::new();
        for (branch_index, branch) in graph.branches().iter().enumerate() {
            thermal_checkpoint(cx)?;
            let pipe = junctions.len() + branch_index;
            for (node, from_node) in [(branch.from, true), (branch.to, false)] {
                if let Some(junction) = slots[node] {
                    connections.push(Connection {
                        from: if from_node { junction } else { pipe },
                        to: if from_node { pipe } else { junction },
                        branch: branch_index,
                        sign: 1.0,
                    });
                } else {
                    let first = incidents[node][0];
                    if first.branch == branch_index {
                        connections.push(Connection {
                            from: pipe,
                            to: junctions.len() + incidents[node][1].branch,
                            branch: first.branch,
                            sign: first.inward_sign,
                        });
                    }
                }
            }
        }
        Ok(Self {
            nodes: nodes.to_vec(),
            junctions,
            connections,
            splices,
        })
    }

    /// A pressure or injection port needs an actual mixing inventory, even if
    /// its current injection happens to vanish. Feedback/scheduled drives use
    /// this same admission on every candidate, before native momentum advances.
    pub(super) fn validate_ports<E>(
        &self,
        cx: &Cx<'_>,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
    ) -> Result<(), ImplicitTransportError<E>> {
        for node in boundaries
            .iter()
            .map(|port| port.node)
            .chain(injections.iter().map(|port| port.node))
        {
            thermal_checkpoint(cx)?;
            if self.nodes.get(node).is_some_and(Option::is_none) {
                return Err(ImplicitTransportError::Invalid(
                    "zero-volume Carreau splice cannot carry pressure or injection ports",
                ));
            }
        }
        Ok(())
    }

    /// Use exactly the original Advected fixed-mass roundoff bound, applied to
    /// the actual rho-scaled incident rates before removing a zero-volume node.
    /// No averaging, projection or hydraulic-tolerance repair is performed.
    /// Each direct link carries the first incident branch's original rate with
    /// its declared orientation. The native binding independently rechecks all
    /// remaining physical cells, so cumulative mismatch cannot become leakage.
    pub(super) fn validate_flows<E>(
        &self,
        cx: &Cx<'_>,
        density: Density,
        volume_flows: &[f64],
    ) -> Result<(), ImplicitTransportError<E>> {
        for splice in &self.splices {
            thermal_checkpoint(cx)?;
            let mut incoming = 0.0;
            let mut outgoing = 0.0;
            for incident in splice.incidents {
                let q = volume_flows[incident.branch];
                let mass = density.value() * q;
                if !mass.is_finite() || (q != 0.0 && mass == 0.0) {
                    return Err(ImplicitTransportError::Invalid(
                        "unrepresentable Carreau splice mass rate",
                    ));
                }
                let inward = incident.inward_sign * mass;
                incoming += inward.max(0.0);
                outgoing += (-inward).max(0.0);
            }
            if !incoming.is_finite()
                || !outgoing.is_finite()
                || (incoming - outgoing).abs() > 64.0 * f64::EPSILON * incoming.max(outgoing)
            {
                return Err(ImplicitTransportError::Invalid(
                    "zero-volume Carreau splice requires native fixed-mass flow continuity",
                ));
            }
        }
        Ok(())
    }
}
