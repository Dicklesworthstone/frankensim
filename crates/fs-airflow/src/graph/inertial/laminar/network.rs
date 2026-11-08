//! Unsteady incompressible networks of the ORIGINAL developed laminar pipes.
//!
//! A boundary-rooted circulation basis enforces every internal-node continuity
//! equation at the level of the continuous dynamics. Only duplicate mean-flow
//! coordinates are eliminated. Every original radial-memory arm, material
//! source and positive viscous source square remains present.
//!
//! With x_mean = U*y and U^T U = I, the native energy is H=||x||^2/2,
//! R_network = T^T R_pipes T and G_network = T^T G_pipes. Reservoir pressures
//! and explicit branch pressure heads have their ORIGINAL flow-conjugate ports.
//! Native exact held-input and continuous-tone owners advance this coupled system.
//! Incidence interconnection: van der Schaft and Maschke (2013),
//! <https://epubs.siam.org/doi/10.1137/110840091>.
//!
//! This is a bounded developed-flow network, not spatial CFD. Junctions have
//! zero storage and add no entrance, mixing or area-change losses. Reference
//! density is common and fixed; every pipe retains its original material,
//! sampled frequency band and Reynolds ceiling. A pressure node is required
//! in each component; one node can serve only as the gauge of a closed loop.
//! Pressure heads are supplied physical sources, not automatically inferred
//! from a geometry or a material. Their work is separate from reservoir work.

mod continuity;
mod prepared;
mod tones;
pub use prepared::{
    LaminarNetworkFlowControlStep, LaminarNetworkFlowPort, LaminarNetworkFlowTolerance,
    PreparedLaminarNetwork, PreparedLaminarNetworkFlowControl,
};
pub use tones::PreparedLaminarNetworkTones;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fluid_transport::FluidTransport;
use fs_math::det;
use fs_phs::{PortHamiltonian, PreparedStepRecord, QuadraticStorage, ResistanceWitness};
use fs_qty::Pressure;
use fs_sparse::Coo;
use fs_sparse::direct::{DirectOrdering, LdltOptions, SymbolicLdlt};

use super::{LaminarPipe, LaminarPipeError, LaminarPipeState, checkpoint, energy, positive};
use continuity::{ContinuityBasis, ROUNDOFF, norm};

/// The existing exact owner reserves 2*(states+ports+1) coordinates, capped at
/// 256. Branch pressure-head ports are included, never hidden from the budget.
pub const MAX_LAMINAR_NETWORK_AUGMENTED: usize = 128;

/// Authored topology. Branch order is the order of the original pipe owners.
/// Pressure-node order is the order of reservoir inputs and boundary outputs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LaminarNetworkTopology {
    pub node_count: usize,
    pub branch_nodes: Vec<[usize; 2]>,
    pub pressure_nodes: Vec<usize>,
}

/// Original energy coordinates, with all internal continuity already eliminated.
/// No pressure multipliers or fictitious pressure energy enter this history.
#[derive(Clone, Debug)]
pub struct LaminarNetworkState {
    coordinates: Vec<f64>,
    accepted_time_s: f64,
    identity: Arc<()>,
}
impl LaminarNetworkState {
    pub const fn accepted_time_s(&self) -> f64 {
        self.accepted_time_s
    }
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        energy(&self.coordinates)
    }
}

/// Actual branch observations of the single coupled trajectory.
#[derive(Clone, Copy, Debug, Default)]
pub struct LaminarNetworkBranchStep {
    pub initial_volume_flow_m3_s: f64,
    pub final_volume_flow_m3_s: f64,
    pub average_volume_flow_m3_s: f64,
    pub displaced_volume_m3: f64,
    pub initial_mean_kinetic_energy_j: f64,
    pub final_mean_kinetic_energy_j: f64,
    pub initial_profile_memory_energy_j: f64,
    pub final_profile_memory_energy_j: f64,
    /// Integral of this pipe's original positive viscous source squares.
    pub viscous_heat_j: f64,
    /// Conservative whole-interval ceiling, not an observed transition.
    pub reynolds_envelope: f64,
}

/// Whole-network native work, storage and original local heat.
#[derive(Clone, Copy, Debug)]
pub struct LaminarNetworkStep {
    pub duration_s: f64,
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    /// Work from reservoir pressures relative to each component's own gauge.
    pub reservoir_work_j: f64,
    /// Work from the explicit pressure heads along authored branch directions.
    pub branch_head_work_j: f64,
    pub maximum_internal_imbalance_m3_s: f64,
    pub maximum_component_boundary_imbalance_m3_s: f64,
    /// The original exact owner integrates total work and heat independently.
    pub native: PreparedStepRecord,
}

/// Source-bearing, constrained native pipe dynamics. Cold construction prepares
/// geometry and steady response; fixed-clock preparation owns trajectory caches.
pub struct LaminarNetwork<M> {
    topology: LaminarNetworkTopology,
    branches: Vec<LaminarPipe<M>>,
    basis: ContinuityBasis,
    offsets: Vec<usize>,
    flow_rows: Vec<Vec<f64>>,
    flow_norms: Vec<f64>,
    source_rows: Vec<Vec<f64>>,
    steady_response: Vec<f64>,
    system: PortHamiltonian,
    identity: Arc<()>,
}

impl<M: FluidTransport> LaminarNetwork<M> {
    pub fn try_new(
        cx: &Cx<'_>,
        topology: LaminarNetworkTopology,
        branches: Vec<LaminarPipe<M>>,
    ) -> Result<Self, LaminarPipeError> {
        checkpoint(cx)?;
        let Some(first) = branches.first() else {
            return Err(LaminarPipeError::Invalid(
                "a laminar network needs original pipe branches",
            ));
        };
        if topology.branch_nodes.len() != branches.len() || branches.len() > 128 {
            return Err(LaminarPipeError::Invalid(
                "one original pipe per bounded authored branch required",
            ));
        }
        let density = first.transport.density_kg_m3();
        let mut arms = 0usize;
        for pipe in &branches {
            checkpoint(cx)?;
            if pipe.transport.density_kg_m3() != density {
                return Err(LaminarPipeError::Invalid(
                    "network branches require exactly common reference density",
                ));
            }
            arms = arms
                .checked_add(pipe.memory.series_terms().len())
                .ok_or(LaminarPipeError::Invalid("network profile extent overflow"))?;
            if arms > MAX_LAMINAR_NETWORK_AUGMENTED {
                return Err(LaminarPipeError::Invalid(
                    "complete network viscous memory exceeds native bound",
                ));
            }
        }
        let scales: Vec<_> = branches.iter().map(|p| p.flow_scales[0]).collect();
        let basis = continuity::compile(
            topology.node_count,
            &topology.branch_nodes,
            &topology.pressure_nodes,
            &scales,
        )
        .map_err(LaminarPipeError::Invalid)?;
        let n = basis.columns + arms;
        let ports = topology.pressure_nodes.len() + branches.len();
        if n == 0 || n + ports + 1 > MAX_LAMINAR_NETWORK_AUGMENTED {
            return Err(LaminarPipeError::Invalid(
                "complete network states and physical pressure ports exceed native exponential bound",
            ));
        }
        let mut flow_rows = Vec::with_capacity(branches.len());
        let mut flow_norms = Vec::with_capacity(branches.len());
        let mut source_rows = Vec::with_capacity(branches.len());
        let mut offsets = Vec::with_capacity(branches.len() + 1);
        let mut resistance = vec![0.0; n * n];
        let mut all_sources = Vec::new();
        let mut input = vec![0.0; n * ports];
        let mut offset = basis.columns;
        for (edge, pipe) in branches.iter().enumerate() {
            checkpoint(cx)?;
            offsets.push(offset);
            let mean = &basis.energy_rows[edge * basis.columns..(edge + 1) * basis.columns];
            let mut flow = vec![0.0; n];
            for (q, u) in flow.iter_mut().zip(mean) {
                *q = u / scales[edge];
            }
            let flow_norm = norm(&flow);
            if !flow_norm.is_finite() {
                return Err(LaminarPipeError::Invalid(
                    "network branch flow map is unrepresentable",
                ));
            }
            let [from, to] = topology.branch_nodes[edge];
            for i in 0..n {
                if let Some(port) = basis.boundary_slot[from] {
                    input[i * ports + port] += flow[i];
                }
                if let Some(port) = basis.boundary_slot[to] {
                    input[i * ports + port] -= flow[i];
                }
                // A REAL pressure-head source port also observes the original
                // branch flow. Unused heads are zero; no observer exponentials.
                input[i * ports + topology.pressure_nodes.len() + edge] = flow[i];
            }
            let mut local = Vec::new();
            let mut row = vec![0.0; n];
            let drag = det::sqrt(pipe.memory.series_drag_per_s());
            for (r, u) in row.iter_mut().zip(mean) {
                *r = drag * u;
            }
            if row.iter().any(|v| *v != 0.0) {
                local.extend_from_slice(&row);
            }
            for (j, term) in pipe.memory.series_terms().iter().enumerate() {
                row.fill(0.0);
                let a = det::sqrt(term.gain_per_s);
                let b = det::sqrt(term.rate_per_s);
                for (r, u) in row.iter_mut().zip(mean) {
                    *r = a * u;
                }
                row[offset + j] = -b;
                local.extend_from_slice(&row);
            }
            if local.len() > n * n || local.iter().any(|v| !v.is_finite()) {
                return Err(LaminarPipeError::Invalid(
                    "complete branch viscous source squares exceed native bound",
                ));
            }
            for row in local.chunks_exact(n) {
                for i in 0..n {
                    for j in 0..=i {
                        let value = row[i] * row[j];
                        resistance[i * n + j] += value;
                        if i != j {
                            resistance[j * n + i] += value;
                        }
                    }
                }
            }
            all_sources.extend_from_slice(&local);
            source_rows.push(local);
            flow_rows.push(flow);
            flow_norms.push(flow_norm);
            offset += pipe.memory.series_terms().len();
        }
        offsets.push(offset);
        let witness = ResistanceWitness::try_from_source_rows(n, resistance.clone(), all_sources)?;
        let mut storage = vec![0.0; n * n];
        for i in 0..n {
            storage[i * n + i] = 1.0;
        }
        let steady_response = steady_response(cx, n, ports, &resistance, &input)?;
        let system = PortHamiltonian::new(
            n,
            ports,
            vec![0.0; n * n],
            resistance,
            input,
            Box::new(QuadraticStorage::new(storage, n)?),
        )?
        .with_resistance_witness(witness)?;
        checkpoint(cx)?;
        Ok(Self {
            topology,
            branches,
            basis,
            offsets,
            flow_rows,
            flow_norms,
            source_rows,
            steady_response,
            system,
            identity: Arc::new(()),
        })
    }

    pub fn topology(&self) -> &LaminarNetworkTopology {
        &self.topology
    }
    pub fn branches(&self) -> &[LaminarPipe<M>] {
        &self.branches
    }
    /// Reservoir ports first, then one signed pressure-head/flow port per pipe.
    /// External composition owns the same source/regime and gauge admission.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn independent_mean_flows(&self) -> usize {
        self.basis.columns
    }

    /// Fixed-density material continuation retains every original pipe identity
    /// and radial history. It never reinitializes the coupled accepted state.
    pub fn with_sources<N: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        sources: Vec<N>,
    ) -> Result<LaminarNetwork<N>, LaminarPipeError> {
        checkpoint(cx)?;
        if sources.len() != self.branches.len() {
            return Err(LaminarPipeError::Invalid(
                "one source per original network branch required",
            ));
        }
        let branches = self
            .branches
            .iter()
            .zip(sources)
            .map(|(pipe, source)| pipe.with_source(cx, source))
            .collect::<Result<Vec<_>, _>>()?;
        let mut next = LaminarNetwork::try_new(cx, self.topology.clone(), branches)?;
        if next.basis.energy_rows != self.basis.energy_rows || next.offsets != self.offsets {
            return Err(LaminarPipeError::Invalid(
                "material continuation changed original circulation coordinates",
            ));
        }
        next.identity = Arc::clone(&self.identity);
        Ok(next)
    }

    pub fn rest_state(
        &self,
        accepted_time_s: f64,
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        self.state(vec![0.0; self.system.state_dim()], accepted_time_s)
    }

    /// Original steady solution for the held reservoir pressures and zero heads.
    pub fn steady_state(
        &self,
        pressures: &[Pressure],
        accepted_time_s: f64,
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        self.steady_state_impl(pressures, None, accepted_time_s)
    }

    pub fn steady_state_with_heads(
        &self,
        pressures: &[Pressure],
        heads: &[Pressure],
        accepted_time_s: f64,
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        self.steady_state_impl(pressures, Some(heads), accepted_time_s)
    }

    fn steady_state_impl(
        &self,
        pressures: &[Pressure],
        heads: Option<&[Pressure]>,
        time: f64,
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        let ports = self.port_count();
        let mut input = vec![0.0; ports];
        self.fill_inputs(pressures, heads, &mut input)?;
        let coordinates: Vec<_> = self
            .steady_response
            .chunks_exact(ports)
            .map(|row| dot(row, &input))
            .collect();
        let state = self.state(coordinates, time)?;
        for edge in 0..self.branches.len() {
            self.branches[edge].admit_flow(dot(&self.flow_rows[edge], &state.coordinates).abs())?;
        }
        Ok(state)
    }

    /// Import original compatible pipe states. Internal volume continuity and
    /// mean-state projection must both fit the explicit ROUNDOFF*branch_count
    /// relative roundoff budget. A physically incompatible state is refused,
    /// never repaired. Every memory coordinate is copied exactly.
    pub fn state_from_branches(
        &self,
        branches: &[LaminarPipeState],
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        if branches.len() != self.branches.len() {
            return Err(LaminarPipeError::Invalid(
                "one original history per network branch required",
            ));
        }
        let time = branches[0].accepted_time_s;
        let mut mean = Vec::with_capacity(branches.len());
        let mut outflow = vec![0.0; self.topology.node_count];
        let mut scale = vec![0.0; self.topology.node_count];
        for (edge, (state, pipe)) in branches.iter().zip(&self.branches).enumerate() {
            if !Arc::ptr_eq(&state.identity, &pipe.identity)
                || state.accepted_time_s != time
                || state.coordinates.len() != pipe.flow_scales.len()
            {
                return Err(LaminarPipeError::Invalid(
                    "original branch histories at one common clock required",
                ));
            }
            let q = state.volume_flow_m3_s();
            pipe.admit_flow(q.abs())?;
            mean.push(state.coordinates[0]);
            let [from, to] = self.topology.branch_nodes[edge];
            outflow[from] += q;
            outflow[to] -= q;
            scale[from] += q.abs();
            scale[to] += q.abs();
        }
        let tolerance = ROUNDOFF * self.branches.len() as f64;
        for node in 0..self.topology.node_count {
            if self.basis.boundary_slot[node].is_none()
                && outflow[node].abs() > tolerance * scale[node]
            {
                return Err(LaminarPipeError::Invalid(
                    "imported branch flows violate internal volume continuity",
                ));
            }
        }
        let r = self.basis.columns;
        let mut coordinates = vec![0.0; self.system.state_dim()];
        for j in 0..r {
            coordinates[j] = (0..branches.len())
                .map(|e| self.basis.energy_rows[e * r + j] * mean[e])
                .sum();
        }
        let difference: Vec<_> = (0..branches.len())
            .map(|e| {
                dot(
                    &self.basis.energy_rows[e * r..(e + 1) * r],
                    &coordinates[..r],
                ) - mean[e]
            })
            .collect();
        if norm(&difference) > tolerance * norm(&mean) {
            return Err(LaminarPipeError::Invalid(
                "import would change physical mean momentum beyond roundoff",
            ));
        }
        for (edge, branch) in branches.iter().enumerate() {
            coordinates[self.offsets[edge]..self.offsets[edge + 1]]
                .copy_from_slice(&branch.coordinates[1..]);
        }
        self.state(coordinates, time)
    }

    /// Recover the same original native pipe state for a separate physical
    /// consumer. This cold extraction allocates; prepared advancement does not.
    pub fn branch_state(
        &self,
        state: &LaminarNetworkState,
        branch: usize,
    ) -> Result<LaminarPipeState, LaminarPipeError> {
        self.admit_owner(state)?;
        let pipe = self
            .branches
            .get(branch)
            .ok_or(LaminarPipeError::Invalid("network branch index"))?;
        let mut coordinates = Vec::with_capacity(pipe.flow_scales.len());
        coordinates.push(dot(
            &self.basis.energy_rows[branch * self.basis.columns..(branch + 1) * self.basis.columns],
            &state.coordinates[..self.basis.columns],
        ));
        coordinates
            .extend_from_slice(&state.coordinates[self.offsets[branch]..self.offsets[branch + 1]]);
        Ok(LaminarPipeState {
            coordinates,
            mean_flow_scale: pipe.flow_scales[0],
            accepted_time_s: state.accepted_time_s,
            identity: Arc::clone(&pipe.identity),
        })
    }

    pub fn branch_volume_flow_m3_s(
        &self,
        state: &LaminarNetworkState,
        branch: usize,
    ) -> Result<f64, LaminarPipeError> {
        self.admit_owner(state)?;
        let row = self
            .flow_rows
            .get(branch)
            .ok_or(LaminarPipeError::Invalid("network branch index"))?;
        Ok(dot(row, &state.coordinates))
    }

    fn state(
        &self,
        coordinates: Vec<f64>,
        time: f64,
    ) -> Result<LaminarNetworkState, LaminarPipeError> {
        if !time.is_finite()
            || time < 0.0
            || coordinates.len() != self.system.state_dim()
            || coordinates.iter().any(|v| !v.is_finite())
            || !energy(&coordinates).is_finite()
        {
            return Err(LaminarPipeError::Invalid(
                "finite complete network history and clock required",
            ));
        }
        Ok(LaminarNetworkState {
            coordinates,
            accepted_time_s: time,
            identity: Arc::clone(&self.identity),
        })
    }
    fn admit_owner(&self, state: &LaminarNetworkState) -> Result<(), LaminarPipeError> {
        if !Arc::ptr_eq(&state.identity, &self.identity) {
            return Err(LaminarPipeError::Invalid(
                "history belongs to a different laminar network",
            ));
        }
        Ok(())
    }
    fn port_count(&self) -> usize {
        self.topology.pressure_nodes.len() + self.branches.len()
    }
    fn fill_inputs(
        &self,
        pressures: &[Pressure],
        heads: Option<&[Pressure]>,
        input: &mut [f64],
    ) -> Result<(), LaminarPipeError> {
        let boundary = self.topology.pressure_nodes.len();
        if pressures.len() != boundary
            || input.len() != self.port_count()
            || heads.is_some_and(|h| h.len() != self.branches.len())
            || pressures.iter().any(|p| !p.value().is_finite())
            || heads.is_some_and(|h| h.iter().any(|p| !p.value().is_finite()))
        {
            return Err(LaminarPipeError::Invalid(
                "complete finite reservoir pressures and branch heads required",
            ));
        }
        // Gauge BEFORE native matrix products. Equal absolute pressures cannot
        // drive a numerical flow through a projected near-zero incidence sum.
        for port in 0..boundary {
            input[port] =
                pressures[port].value() - pressures[self.basis.pressure_reference[port]].value();
        }
        for edge in 0..self.branches.len() {
            input[boundary + edge] = heads.map_or(0.0, |h| h[edge].value());
        }
        if input.iter().any(|v| !v.is_finite()) {
            return Err(LaminarPipeError::Invalid(
                "pressure differences are unrepresentable",
            ));
        }
        Ok(())
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Reuse the existing sparse direct owner for the cold SPD response. Its
/// positive inertia and the ORIGINAL componentwise residual are both checked.
fn steady_response(
    cx: &Cx<'_>,
    n: usize,
    ports: usize,
    resistance: &[f64],
    input: &[f64],
) -> Result<Vec<f64>, LaminarPipeError> {
    let mut matrix = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            let value = resistance[i * n + j];
            if value != 0.0 {
                matrix.push(i, j, value);
            }
        }
    }
    let matrix = matrix.assemble();
    let symbolic = SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd)
        .map_err(|_| LaminarPipeError::Invalid("network steady resistance pattern refused"))?;
    let factor = symbolic
        .factor(&matrix, &LdltOptions::default())
        .map_err(|_| {
            LaminarPipeError::Invalid("network steady resistance factorization refused")
        })?;
    if factor.inertia().positive != n || factor.inertia().negative != 0 {
        return Err(LaminarPipeError::Invalid(
            "network original viscous resistance is not positive definite",
        ));
    }
    let mut map = vec![0.0; n * ports];
    for port in 0..ports {
        checkpoint(cx)?;
        let rhs: Vec<_> = (0..n).map(|i| input[i * ports + port]).collect();
        let solution = factor.solve(&rhs);
        for i in 0..n {
            let mut value = 0.0;
            let mut scale = rhs[i].abs();
            for j in 0..n {
                let term = resistance[i * n + j] * solution[j];
                value += term;
                scale += term.abs();
            }
            if !solution[i].is_finite()
                || !scale.is_finite()
                || (value - rhs[i]).abs() > ROUNDOFF * n as f64 * scale
            {
                return Err(LaminarPipeError::Invalid(
                    "network steady response does not satisfy original momentum",
                ));
            }
            map[i * ports + port] = solution[i];
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests;
