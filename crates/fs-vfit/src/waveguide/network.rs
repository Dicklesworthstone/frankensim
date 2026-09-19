//! Connected characteristic sections with lossless pressure junctions.
//!
//! At a junction, p = 2 sum(a_i / Z_i) / sum(1 / Z_i), b_i = p - a_i.
//! Thus pressure is continuous, sum((a_i-b_i)/Z_i) = 0, and incident and
//! departing wave powers agree in exact arithmetic. See J. O. Smith,
//! Physical Audio Signal Processing, "Adaptors for Wave Digital Elements",
//! general parallel adaptor. All segment delays are positive: scattering
//! reads the SAME old wave state before any segment shifts, including cycles.
//!
//! This is a fixed graph of lossless uniform sections and ideal zero-volume
//! junctions, with one degree-one inlet and passive resistive, R-L-C or relaxation
//! terminal loads. Reactive boundary storage participates in the total balance.
//! No extra sample of junction delay, fitted-filter energy claim, branch end
//! correction, distributed loss, fractional delay or radiation law is implied.

use super::{PassiveWaveguide, WaveguideError, WaveguideSpec};
use core::mem::size_of;
use crate::impedance::{ImpedanceState, SeriesImpedanceSpec};
use crate::relaxation::RelaxationImpedanceSpec;
mod load;
use load::{TerminalFrame, TerminalLoad};

/// Boundary or ideal connecting junction at a graph node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NetworkNode {
    /// Exactly one node is the externally driven, degree-one inlet.
    Inlet,
    /// At least two incident sections; no added storage or source.
    Junction,
    /// Degree-one memoryless load; pressure reflectance must lie in [-1, 1].
    Termination { reflection: f64 },
    /// Degree-one passive reactive load, initially at zero energy. Its physical
    /// state is retained; memoryless reflection controls cannot erase it.
    Impedance { load: SeriesImpedanceSpec },
    /// Degree-one frequency-dependent boundary loss. Each relaxation branch
    /// retains its own inertive history and participates in the energy balance.
    Relaxation { load: RelaxationImpedanceSpec },
}

/// One lossless section. Endpoint order sets direction, not physical authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NetworkSegment {
    /// Distinct, zero-based node indices. Parallel sections and cycles are allowed.
    pub nodes: [usize; 2],
    /// Positive one-way transit in samples.
    pub one_way_samples: usize,
    /// Positive pressure/volume-flow impedance [Pa s/m^3].
    pub impedance_pa_s_m3: f64,
}

/// Observations at a node on the last ACCEPTED step (initially zero).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct NodeFrame {
    /// Common pressure [Pa].
    pub pressure_pa: f64,
    /// Sum of flows from adjacent sections INTO this node [m^3/s].
    /// Approximately zero for a junction; negative of inlet flow at the inlet.
    pub net_flow_into_node_m3_s: f64,
    /// Irreversibly dissipated terminal energy [J]; zero elsewhere. Reactive
    /// port work may be negative, but is never confused with dissipation.
    pub absorbed_energy_j: f64,
    /// End-of-step energy retained in this terminal's inertance/compliance [J].
    pub stored_energy_j: f64,
    /// Change in terminal storage [J], independently evaluated from its state.
    pub storage_change_j: f64,
}

/// Same-step inlet work, all-terminal dissipation and wave PLUS reactive storage.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct NetworkFrame {
    /// Returning wave used by this step's inlet solve [Pa].
    pub incoming_pressure_pa: f64,
    /// Inlet pressure [Pa].
    pub inlet_pressure_pa: f64,
    /// Volume flow INTO the network [m^3/s].
    pub inlet_flow_m3_s: f64,
    /// End-of-step wave plus reactive boundary energy [J].
    pub stored_energy_j: f64,
    /// Energy in bidirectional delay buffers alone [J].
    pub wave_stored_energy_j: f64,
    /// Energy in all reactive terminal states alone [J].
    pub load_stored_energy_j: f64,
    /// Difference in independently evaluated combined storage [J].
    pub storage_change_j: f64,
    /// Inlet pressure times inlet flow times dt [J].
    pub inlet_work_j: f64,
    /// Sum of resistive terminal dissipations [J], excluding reversible storage.
    pub terminal_loss_j: f64,
    /// Incident minus departing energy at internal junctions [J]. Rounding
    /// diagnostic only: it is NOT treated as dissipation or used to repair storage.
    pub junction_residual_j: f64,
}
impl NetworkFrame {
    /// Total balance [J], without subtracting away any scattering defect.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.terminal_loss_j - self.inlet_work_j
    }
}

/// A fixed topology and retained traveling waves. No allocation during stepping.
/// Work is O(number of ports + sum(log(section delay))) per preview/commit.
/// All reductions follow declared node/section order, independent of callbacks.
pub struct WaveguideNetwork {
    nodes: Vec<NetworkNode>,
    offsets: Vec<usize>,
    ports: Vec<usize>,
    weights: Vec<f64>,
    lines: Vec<PassiveWaveguide>,
    departing: Vec<f64>,
    observed: Vec<NodeFrame>,
    candidate: Vec<NodeFrame>,
    loads: Vec<Option<TerminalLoad>>,
    load_candidate: Vec<TerminalFrame>,
    inlet_port: usize,
    time_step_s: f64,
}

fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, WaveguideError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| WaveguideError("network allocation failed"))?;
    out.resize(n, value);
    Ok(out)
}

fn representative(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    node
}

impl WaveguideNetwork {
    /// Requested heap-payload admission, including record arrays, scratch,
    /// energy trees, and a temporary constructor cursor. Excludes caller input
    /// vectors, allocator overhead and RSS. Arithmetic overflow refuses.
    ///
    /// # Errors
    /// Empty graph, zero/unrepresentable delay or size arithmetic overflow.
    pub fn required_memory_bytes(
        node_count: usize, segments: &[NetworkSegment],
    ) -> Result<usize, WaveguideError> {
        let overflow = WaveguideError("network size overflow");
        if node_count == 0 || segments.is_empty() {
            return Err(WaveguideError("network must contain nodes and sections"));
        }
        let mut bytes = size_of::<usize>(); // offsets has node_count + 1 entries
        for (count, width) in [
            (node_count, size_of::<NetworkNode>() + 2 * size_of::<usize>()
                + 2 * size_of::<NodeFrame>() + size_of::<Option<TerminalLoad>>()
                + size_of::<TerminalFrame>()),
            (segments.len(), size_of::<PassiveWaveguide>() + 2 * size_of::<usize>()
                + 4 * size_of::<f64>()),
        ] {
            bytes = bytes.checked_add(count.checked_mul(width).ok_or(overflow)?).ok_or(overflow)?;
        }
        for segment in segments {
            if segment.one_way_samples == 0 {
                return Err(WaveguideError("every network section needs a positive transit"));
            }
            let leaves = segment.one_way_samples.checked_next_power_of_two().ok_or(overflow)?;
            let values = segment.one_way_samples.checked_mul(2)
                .and_then(|n| leaves.checked_mul(4).and_then(|e| n.checked_add(e)))
                .ok_or(overflow)?;
            bytes = bytes.checked_add(values.checked_mul(size_of::<f64>()).ok_or(overflow)?)
                .ok_or(overflow)?;
        }
        Ok(bytes)
    }

    /// Admit a connected graph, then allocate zero-state sections and scratch.
    /// Inlet/terminal nodes have degree one; internal junctions have degree >=2.
    ///
    /// # Errors
    /// Topology, passive-load, finite impedance/time, memory or allocation refusal.
    pub fn new(
        nodes: &[NetworkNode], segments: &[NetworkSegment], time_step_s: f64,
        max_memory_bytes: usize,
    ) -> Result<Self, WaveguideError> {
        let required = Self::required_memory_bytes(nodes.len(), segments)?;
        if required > max_memory_bytes {
            return Err(WaveguideError("network payload exceeds memory budget"));
        }
        if !time_step_s.is_finite() || time_step_s <= 0.0 {
            return Err(WaveguideError("network time step must be positive and finite"));
        }
        for segment in segments {
            if segment.nodes[0] == segment.nodes[1]
                || segment.nodes.iter().any(|&n| n >= nodes.len())
                || !segment.impedance_pa_s_m3.is_finite() || segment.impedance_pa_s_m3 <= 0.0
            {
                return Err(WaveguideError("network section needs valid distinct nodes and positive finite impedance"));
            }
        }
        let port_count = segments.len().checked_mul(2)
            .ok_or(WaveguideError("network port count overflow"))?;
        let mut offsets = buffer(nodes.len() + 1, 0usize)?;
        for segment in segments {
            for &node in &segment.nodes { offsets[node + 1] += 1; }
        }
        let mut inlet_node = None;
        for (n, kind) in nodes.iter().enumerate() {
            let degree = offsets[n + 1];
            match *kind {
                NetworkNode::Inlet if degree == 1 && inlet_node.is_none() => inlet_node = Some(n),
                NetworkNode::Junction if degree >= 2 => {},
                NetworkNode::Impedance { load } if degree == 1 => { load.validate()?; },
                NetworkNode::Relaxation { .. } if degree == 1 => {},
                NetworkNode::Termination { reflection }
                    if degree == 1 && reflection.is_finite() && reflection.abs() <= 1.0 => {},
                _ => return Err(WaveguideError("invalid node degree, repeated inlet or active/nonfinite termination")),
            }
        }
        let inlet_node = inlet_node.ok_or(WaveguideError("network requires one inlet"))?;
        for n in 0..nodes.len() { offsets[n + 1] += offsets[n]; }
        let mut cursor = buffer(nodes.len(), 0usize)?;
        cursor.copy_from_slice(&offsets[..nodes.len()]);
        let mut ports = buffer(port_count, 0usize)?;
        for (edge, segment) in segments.iter().enumerate() {
            for side in 0..2 {
                let node = segment.nodes[side];
                ports[cursor[node]] = 2 * edge + side;
                cursor[node] += 1;
            }
        }
        // Reuse constructor cursor as union-find storage; disconnected pockets
        // cannot silently become unexcited "success" sections.
        for (n, p) in cursor.iter_mut().enumerate() { *p = n; }
        for segment in segments {
            let a = representative(&mut cursor, segment.nodes[0]);
            let b = representative(&mut cursor, segment.nodes[1]);
            cursor[a] = b;
        }
        let root = representative(&mut cursor, inlet_node);
        if (0..nodes.len()).any(|n| representative(&mut cursor, n) != root) {
            return Err(WaveguideError("all network sections must connect to the inlet"));
        }
        let mut weights = buffer(port_count, 0.0)?;
        for n in 0..nodes.len() {
            let adjacent = &ports[offsets[n]..offsets[n + 1]];
            let min_z = adjacent.iter().map(|&p| segments[p / 2].impedance_pa_s_m3)
                .fold(f64::INFINITY, f64::min);
            let mut sum = 0.0;
            for &port in adjacent {
                weights[port] = min_z / segments[port / 2].impedance_pa_s_m3;
                sum += weights[port];
            }
            for &port in adjacent {
                weights[port] /= sum;
                if !weights[port].is_finite() || weights[port] <= 0.0 {
                    return Err(WaveguideError("junction admittance weights are not representable"));
                }
            }
        }
        let mut lines = Vec::new();
        lines.try_reserve_exact(segments.len()).map_err(|_| WaveguideError("network allocation failed"))?;
        for segment in segments {
            lines.push(PassiveWaveguide::new(WaveguideSpec {
                one_way_samples: segment.one_way_samples,
                impedance_pa_s_m3: segment.impedance_pa_s_m3, time_step_s,
                reflection: 0.0, // endpoint waves belong to the graph, not this unused load
                max_memory_bytes,
            })?);
        }
        let mut retained_nodes = buffer(nodes.len(), NetworkNode::Junction)?;
        retained_nodes.copy_from_slice(nodes);
        let inlet_port = ports[offsets[inlet_node]];
        let mut loads = buffer(nodes.len(), None)?;
        for (n, kind) in nodes.iter().enumerate() {
            let z = segments[ports[offsets[n]] / 2].impedance_pa_s_m3;
            loads[n] = TerminalLoad::new(*kind, z, time_step_s)?;
        }
        Ok(Self {
            nodes: retained_nodes, offsets, ports, weights, lines,
            departing: buffer(port_count, 0.0)?, observed: buffer(nodes.len(), NodeFrame::default())?,
            candidate: buffer(nodes.len(), NodeFrame::default())?, inlet_port, time_step_s,
            loads, load_candidate: buffer(nodes.len(), TerminalFrame::default())?,
        })
    }

    fn arriving(&self, port: usize) -> f64 {
        let line = &self.lines[port / 2];
        if port % 2 == 0 { line.incoming_pressure_pa() } else { line.waves[line.head] }
    }

    /// Returning inlet wave before the next sample advances.
    #[must_use]
    pub fn incoming_pressure_pa(&self) -> f64 { self.arriving(self.inlet_port) }

    /// Impedance of the unique inlet section [Pa s/m^3].
    #[must_use]
    pub fn inlet_impedance_pa_s_m3(&self) -> f64 {
        self.lines[self.inlet_port / 2].spec.impedance_pa_s_m3
    }

    /// Actual delay-buffer storage, summed in declared section order [J].
    #[must_use]
    pub fn wave_stored_energy_j(&self) -> f64 {
        self.lines.iter().fold(0.0, |sum, line| sum + line.stored_energy_j())
    }

    /// Actual reactive storage summed in declared node order [J].
    #[must_use]
    pub fn load_stored_energy_j(&self) -> f64 {
        self.loads.iter().fold(0.0, |sum, load| {
            sum + load.as_ref().map_or(0.0, TerminalLoad::stored_energy_j)
        })
    }

    /// Actual wave and load storage [J]. No accumulated work integral.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.wave_stored_energy_j() + self.load_stored_energy_j()
    }

    /// Accepted base R-L-C coordinates; None for a nonreactive or unknown node.
    /// Relaxation histories are available separately through `terminal_relaxation_flows`.
    #[must_use]
    pub fn terminal_state(&self, node: usize) -> Option<ImpedanceState> {
        self.loads.get(node).and_then(Option::as_ref).map(TerminalLoad::state)
    }

    /// Accepted internal relaxation flows [m^3/s]; None for other node kinds.
    /// Preview, failed samples and refused reflection controls preserve these.
    #[must_use]
    pub fn terminal_relaxation_flows(&self, node: usize) -> Option<&[f64]> {
        self.loads.get(node).and_then(Option::as_ref).and_then(TerminalLoad::relaxation_flows)
    }

    /// Last accepted node observation; a preview never replaces this value.
    #[must_use]
    pub fn node_frame(&self, node: usize) -> Option<&NodeFrame> { self.observed.get(node) }

    /// Current boundary/junction declaration.
    #[must_use]
    pub fn node(&self, node: usize) -> Option<NetworkNode> { self.nodes.get(node).copied() }

    /// Change an ideal memoryless termination BETWEEN samples, preserving all
    /// traveling waves. Every value remains passive; no stored load state exists
    /// to reset. Reactive terminals REFUSE this control: it cannot discard their
    /// stored energy. This is not a moving pad, displacement work, or radiation.
    ///
    /// # Errors
    /// Unknown/nonterminal node or nonfinite/active reflectance; no mutation.
    pub fn set_terminal_reflection(&mut self, node: usize, reflection: f64) -> Result<(), WaveguideError> {
        if !reflection.is_finite() || reflection.abs() > 1.0 {
            return Err(WaveguideError("terminal reflectance must be finite and passive"));
        }
        match self.nodes.get_mut(node) {
            Some(NetworkNode::Termination { reflection: current }) => { *current = reflection; Ok(()) },
            _ => Err(WaveguideError("reflection control requires an existing terminal node")),
        }
    }

    /// Prepare simultaneous scattering and prospective storage. Scratch may
    /// change; delay buffers, energy and accepted observations do not. No heap
    /// allocation occurs. A later step recomputes this preview before committing.
    ///
    /// # Errors
    /// Nonfinite input, port quantities, candidate wave storage or work.
    pub fn preview_step(&mut self, outgoing: f64) -> Result<NetworkFrame, WaveguideError> {
        if !outgoing.is_finite() { return Err(WaveguideError("network input must be finite")); }
        let incoming = self.incoming_pressure_pa();
        let z = self.inlet_impedance_pa_s_m3();
        let mut terminal_loss = 0.0;
        let mut junction_residual = 0.0;
        let mut load_stored = 0.0;
        for n in 0..self.nodes.len() {
            let range = self.offsets[n]..self.offsets[n + 1];
            let first = self.ports[range.start];
            let pressure = match self.nodes[n] {
                NetworkNode::Inlet => outgoing + self.arriving(first),
                NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. } => {
                    let a = self.arriving(first);
                    let trial = self.loads[n].as_ref().expect("admitted reactive terminal").preview_step(a)?;
                    self.load_candidate[n] = trial;
                    load_stored += trial.port().stored_energy_j;
                    trial.port().pressure_pa
                },
                NetworkNode::Termination { reflection } => {
                    let a = self.arriving(first);
                    a + reflection * a
                },
                NetworkNode::Junction => {
                    let mut mean = 0.0;
                    for index in range.clone() {
                        let port = self.ports[index];
                        mean += self.weights[port] * self.arriving(port);
                    }
                    2.0 * mean
                },
            };
            let mut observation = NodeFrame { pressure_pa: pressure, ..NodeFrame::default() };
            for index in range {
                let port = self.ports[index];
                let a = self.arriving(port);
                let b = match self.nodes[n] {
                    NetworkNode::Inlet => outgoing,
                    NetworkNode::Termination { reflection } => reflection * a,
                    NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. } => self.load_candidate[n].port().reflected_pressure_pa,
                    NetworkNode::Junction => pressure - a,
                };
                let line = &self.lines[port / 2];
                let flow = (a - b) / line.spec.impedance_pa_s_m3;
                if !b.is_finite() || !flow.is_finite() {
                    return Err(WaveguideError("network scattering left the finite set"));
                }
                observation.net_flow_into_node_m3_s += flow;
                match self.nodes[n] {
                    NetworkNode::Termination { reflection } => {
                        observation.absorbed_energy_j = line.wave_energy(a)
                            * (1.0 - reflection) * (1.0 + reflection);
                        terminal_loss += observation.absorbed_energy_j;
                    },
                    NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. } => {
                        let trial = self.load_candidate[n].port();
                        observation.absorbed_energy_j = trial.dissipated_energy_j;
                        observation.stored_energy_j = trial.stored_energy_j;
                        observation.storage_change_j = trial.storage_change_j;
                        terminal_loss += trial.dissipated_energy_j;
                    },
                    NetworkNode::Junction => junction_residual += line.wave_energy(a) - line.wave_energy(b),
                    NetworkNode::Inlet => {},
                }
                self.departing[port] = b;
            }
            if ![pressure, observation.net_flow_into_node_m3_s, observation.absorbed_energy_j]
                .iter().all(|x| x.is_finite())
            {
                return Err(WaveguideError("network node observation left the finite set"));
            }
            self.candidate[n] = observation;
        }
        let mut stored = 0.0;
        for (i, line) in self.lines.iter().enumerate() {
            stored += line.replaced_root(0, line.wave_energy(self.departing[2 * i]))
                + line.replaced_root(2 * line.leaves, line.wave_energy(self.departing[2 * i + 1]));
        }
        let wave_stored = stored;
        stored += load_stored;
        let frame = NetworkFrame {
            incoming_pressure_pa: incoming, inlet_pressure_pa: outgoing + incoming,
            inlet_flow_m3_s: (outgoing - incoming) / z,
            stored_energy_j: stored, storage_change_j: stored - self.stored_energy_j(),
            wave_stored_energy_j: wave_stored, load_stored_energy_j: load_stored,
            inlet_work_j: (outgoing + incoming) * (((outgoing - incoming) / z) * self.time_step_s),
            terminal_loss_j: terminal_loss, junction_residual_j: junction_residual,
        };
        if ![frame.inlet_pressure_pa, frame.inlet_flow_m3_s, stored, frame.storage_change_j,
            frame.inlet_work_j, terminal_loss, junction_residual, frame.balance_residual_j()]
            .iter().all(|x| x.is_finite())
        {
            return Err(WaveguideError("network storage or work left the finite set"));
        }
        Ok(frame)
    }

    /// Scatter every old endpoint, validate every candidate, THEN advance all
    /// sections. Refusal advances no physical state or accepted observation.
    ///
    /// # Errors
    /// Same finite-set checks as `preview_step`; no allocation in this call.
    pub fn step(&mut self, outgoing: f64) -> Result<NetworkFrame, WaveguideError> {
        let frame = self.preview_step(outgoing)?;
        for (i, line) in self.lines.iter_mut().enumerate() {
            line.commit_pair(self.departing[2 * i], self.departing[2 * i + 1]);
        }
        for (load, candidate) in self.loads.iter_mut().zip(&self.load_candidate) {
            if let Some(load) = load { load.accept_frame(*candidate); }
        }
        self.observed.copy_from_slice(&self.candidate);
        Ok(frame)
    }
}

#[cfg(test)]
mod tests;
