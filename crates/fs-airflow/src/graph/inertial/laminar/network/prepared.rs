//! One native coupled trajectory and original per-pipe source heat.
use super::*;
mod harmonic;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh, PreparedSourceLoss,
};

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Prepare once for a positive held-pressure interval. All branch heat
    /// caches observe the SAME projected generator and original source rows.
    /// Successful warm advancement allocates nothing and solves no new system.
    pub fn prepare(
        &self,
        cx: &Cx<'_>,
        duration_s: f64,
    ) -> Result<PreparedLaminarNetwork<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        if !positive(duration_s) {
            return Err(LaminarPipeError::Invalid(
                "positive finite network interval required",
            ));
        }
        let n = self.system.state_dim();
        let ports = self.port_count();
        let native = PreparedExactZoh::new(
            &self.system,
            duration_s,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: 1.0,
            },
        )?;
        let (_, resistance, input) = self.system.structure();
        let generator: Vec<_> = resistance.iter().map(|v| -v).collect();
        let mut heat = Vec::with_capacity(self.branches.len());
        for source in &self.source_rows {
            checkpoint(cx)?;
            heat.push(
                PreparedSourceLoss::new_controlled(
                    n,
                    ports,
                    &generator,
                    input,
                    source,
                    duration_s,
                    DenseExponentialBudget {
                        max_dimension: 2 * n,
                        max_dense_bytes: usize::MAX,
                        max_multiply_terms: u64::MAX,
                    },
                    || cx.checkpoint().is_err(),
                )
                .map_err(prepared_error)?,
            );
        }
        checkpoint(cx)?;
        Ok(PreparedLaminarNetwork {
            network: self,
            native,
            heat,
            pending: vec![0.0; n],
            input: vec![0.0; ports],
            output: vec![0.0; ports],
            steady: vec![0.0; n],
            drive: vec![0.0; n],
            difference: vec![0.0; n],
            residual: vec![0.0; n],
            heat_scratch: vec![0.0; 2 * n],
            records: vec![LaminarNetworkBranchStep::default(); self.branches.len()],
            pending_records: vec![LaminarNetworkBranchStep::default(); self.branches.len()],
            boundary_flows: vec![0.0; self.topology.pressure_nodes.len()],
            flow_audit: NetworkFlowAudit::new(self),
            last: None,
        })
    }
}

/// Prepared native state/heat transition with reusable transactional scratch.
pub struct PreparedLaminarNetwork<'a, M> {
    network: &'a LaminarNetwork<M>,
    native: PreparedExactZoh<'a>,
    heat: Vec<PreparedSourceLoss>,
    pending: Vec<f64>,
    input: Vec<f64>,
    output: Vec<f64>,
    steady: Vec<f64>,
    drive: Vec<f64>,
    difference: Vec<f64>,
    residual: Vec<f64>,
    heat_scratch: Vec<f64>,
    records: Vec<LaminarNetworkBranchStep>,
    pending_records: Vec<LaminarNetworkBranchStep>,
    boundary_flows: Vec<f64>,
    flow_audit: NetworkFlowAudit,
    last: Option<LaminarNetworkStep>,
}

impl<M: FluidTransport> PreparedLaminarNetwork<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }
    /// Empty before first acceptance; unchanged on any later refusal.
    pub fn branch_steps(&self) -> &[LaminarNetworkBranchStep] {
        if self.last.is_some() {
            &self.records
        } else {
            &[]
        }
    }
    /// Positive flow is from each original reservoir INTO the network.
    /// These are native trajectory averages, never endpoint averages.
    pub fn boundary_volume_flows_m3_s(&self) -> &[f64] {
        if self.last.is_some() {
            &self.boundary_flows
        } else {
            &[]
        }
    }
    pub const fn last_step(&self) -> Option<LaminarNetworkStep> {
        self.last
    }

    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarNetworkState,
        pressures: &[Pressure],
    ) -> Result<LaminarNetworkStep, LaminarPipeError> {
        self.advance_impl(cx, state, pressures, None)
    }
    /// Explicit physical pressure rises along authored branch directions.
    /// Their work is pressure head times the ORIGINAL integrated branch flow.
    pub fn advance_with_heads(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarNetworkState,
        pressures: &[Pressure],
        heads: &[Pressure],
    ) -> Result<LaminarNetworkStep, LaminarPipeError> {
        self.advance_impl(cx, state, pressures, Some(heads))
    }

    fn advance_impl(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarNetworkState,
        pressures: &[Pressure],
        heads: Option<&[Pressure]>,
    ) -> Result<LaminarNetworkStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.network.admit_owner(state)?;
        self.network
            .fill_inputs(pressures, heads, &mut self.input)?;
        let dt = self.duration_s();
        let next_time = state.accepted_time_s + dt;
        if !next_time.is_finite() || next_time <= state.accepted_time_s {
            return Err(LaminarPipeError::Invalid(
                "network physical clock must advance finitely",
            ));
        }
        let n = self.pending.len();
        let ports = self.input.len();
        let boundary = self.network.topology.pressure_nodes.len();
        let (_, resistance, g) = self.network.system.structure();
        for i in 0..n {
            self.steady[i] = dot(
                &self.network.steady_response[i * ports..(i + 1) * ports],
                &self.input,
            );
            let row = &g[i * ports..(i + 1) * ports];
            self.drive[i] = dot(row, &self.input);
            let magnitude: f64 = row
                .iter()
                .zip(&self.input)
                .map(|(g, u)| (g * u).abs())
                .sum();
            self.residual[i] = self.drive[i].abs() + ROUNDOFF * ports as f64 * magnitude;
            self.difference[i] = state.coordinates[i] - self.steady[i];
        }
        let initial_norm = norm(&state.coordinates);
        // Opposing pressure sources can nearly cancel in Gu. Bound their
        // uncancelled dot-product roundoff before taking the forcing norm.
        let transient = initial_norm + dt * norm(&self.residual);
        // Passive contraction about x_infinity gives a useful long-interval
        // bound even at steady flow. Retain the actual solve residual and a
        // roundoff allowance, so a rounded steady response is not treated as
        // an exact equilibrium. Also retain the independent forcing bound.
        for i in 0..n {
            let row = &resistance[i * n..(i + 1) * n];
            let restoring = dot(row, &self.steady);
            let magnitude: f64 = row
                .iter()
                .zip(&self.steady)
                .map(|(r, x)| (r * x).abs())
                .sum();
            let drive_magnitude: f64 = g[i * ports..(i + 1) * ports]
                .iter()
                .zip(&self.input)
                .map(|(g, u)| (g * u).abs())
                .sum();
            self.residual[i] = (self.drive[i] - restoring).abs()
                + ROUNDOFF * (n as f64 * magnitude + ports as f64 * drive_magnitude);
        }
        let deviation = norm(&self.difference) + dt * norm(&self.residual);
        let shifted = norm(&self.steady) + deviation;
        let inflation = 1.0 + ROUNDOFF * n as f64;
        let coordinate_bound = transient.min(shifted) * inflation;
        if !transient.is_finite() || !deviation.is_finite() || !coordinate_bound.is_finite() {
            return Err(LaminarPipeError::Invalid(
                "network whole-interval momentum envelope is unrepresentable",
            ));
        }
        for edge in 0..self.network.branches.len() {
            checkpoint(cx)?;
            let row = &self.network.flow_rows[edge];
            let steady_magnitude: f64 = row
                .iter()
                .zip(&self.steady)
                .map(|(q, x)| (q * x).abs())
                .sum();
            let observed_steady =
                dot(row, &self.steady).abs() + ROUNDOFF * n as f64 * steady_magnitude;
            let row_norm = self.network.flow_norms[edge];
            let bound =
                (row_norm * transient).min(observed_steady + row_norm * deviation) * inflation;
            if !bound.is_finite() {
                return Err(LaminarPipeError::Invalid(
                    "network branch Reynolds envelope is unrepresentable",
                ));
            }
            self.pending_records[edge].reynolds_envelope =
                self.network.branches[edge].admit_flow(bound)?;
        }
        let native = self.native.step_into(
            &state.coordinates,
            &self.input,
            &mut self.pending,
            &mut self.output,
        )?;
        let mut local_heat = 0.0;
        for (edge, heat) in self.heat.iter().enumerate() {
            checkpoint(cx)?;
            let dissipated = heat
                .evaluate_into_controlled(
                    &state.coordinates,
                    &self.input,
                    &mut self.heat_scratch,
                    || cx.checkpoint().is_err(),
                )
                .map_err(prepared_error)?;
            local_heat += dissipated;
            let row = &self.network.flow_rows[edge];
            let r = self.network.basis.columns;
            let mean = &self.network.basis.energy_rows[edge * r..(edge + 1) * r];
            let initial_mean = dot(mean, &state.coordinates[..r]);
            let final_mean = dot(mean, &self.pending[..r]);
            let start = self.network.offsets[edge];
            let end = self.network.offsets[edge + 1];
            let report = &mut self.pending_records[edge];
            report.initial_volume_flow_m3_s = dot(row, &state.coordinates);
            report.final_volume_flow_m3_s = dot(row, &self.pending);
            report.average_volume_flow_m3_s = self.output[boundary + edge];
            report.displaced_volume_m3 = dt * report.average_volume_flow_m3_s;
            report.initial_mean_kinetic_energy_j = 0.5 * initial_mean * initial_mean;
            report.final_mean_kinetic_energy_j = 0.5 * final_mean * final_mean;
            report.initial_profile_memory_energy_j = energy(&state.coordinates[start..end]);
            report.final_profile_memory_energy_j = energy(&self.pending[start..end]);
            report.viscous_heat_j = dissipated;
            if ![
                report.initial_volume_flow_m3_s,
                report.final_volume_flow_m3_s,
                report.average_volume_flow_m3_s,
                report.displaced_volume_m3,
                report.initial_mean_kinetic_energy_j,
                report.final_mean_kinetic_energy_j,
                report.initial_profile_memory_energy_j,
                report.final_profile_memory_energy_j,
                report.viscous_heat_j,
            ]
            .into_iter()
            .all(f64::is_finite)
                || report.viscous_heat_j < 0.0
                || (report.average_volume_flow_m3_s != 0.0 && report.displaced_volume_m3 == 0.0)
            {
                return Err(LaminarPipeError::Invalid(
                    "network branch flow, heat or storage is unrepresentable",
                ));
            }
        }
        let heat_allowance = ROUNDOFF * n as f64 * local_heat.abs().max(native.dissipated.abs());
        if !local_heat.is_finite() || (local_heat - native.dissipated).abs() > heat_allowance {
            return Err(LaminarPipeError::Invalid(
                "original branch heat does not close the native network account",
            ));
        }
        let endpoint = self.audit_flows(false, coordinate_bound)?;
        let average = self.audit_flows(true, coordinate_bound)?;
        let reservoir_work_j: f64 = (0..boundary)
            .map(|p| self.input[p] * self.output[p] * dt)
            .sum();
        let branch_head_work_j: f64 = (boundary..ports)
            .map(|p| self.input[p] * self.output[p] * dt)
            .sum();
        // Recovery and supply can cancel inside either port group. Admission
        // uses every uncancelled work term; signed reports remain unchanged.
        let work_scale = (0..ports)
            .map(|p| (self.input[p] * self.output[p] * dt).abs())
            .sum::<f64>()
            + native.supplied.abs();
        let report = LaminarNetworkStep {
            duration_s: dt,
            initial_kinetic_energy_j: state.kinetic_energy_j(),
            final_kinetic_energy_j: energy(&self.pending),
            reservoir_work_j,
            branch_head_work_j,
            maximum_internal_imbalance_m3_s: endpoint.0.max(average.0),
            maximum_component_boundary_imbalance_m3_s: endpoint.1.max(average.1),
            native,
        };
        if ![
            report.initial_kinetic_energy_j,
            report.final_kinetic_energy_j,
            reservoir_work_j,
            branch_head_work_j,
            work_scale,
        ]
        .into_iter()
        .all(f64::is_finite)
            || (reservoir_work_j + branch_head_work_j - native.supplied).abs()
                > ROUNDOFF * n as f64 * work_scale
        {
            return Err(LaminarPipeError::Invalid(
                "original reservoir and head work do not close native supply",
            ));
        }
        checkpoint(cx)?;
        state.coordinates.copy_from_slice(&self.pending);
        state.accepted_time_s = next_time;
        core::mem::swap(&mut self.records, &mut self.pending_records);
        self.boundary_flows
            .copy_from_slice(&self.output[..boundary]);
        self.last = Some(report);
        Ok(report)
    }

    fn audit_flows(
        &mut self,
        average: bool,
        coordinate_bound: f64,
    ) -> Result<(f64, f64), LaminarPipeError> {
        self.flow_audit.audit(
            self.network,
            &self.pending_records,
            &self.output,
            average,
            coordinate_bound,
        )
    }
}

/// The same original topology audit serves held and continuous trajectories.
pub(super) struct NetworkFlowAudit {
    node_flow: Vec<f64>,
    node_scale: Vec<f64>,
    component_flow: Vec<f64>,
    component_scale: Vec<f64>,
}

impl NetworkFlowAudit {
    pub(super) fn new<M: FluidTransport>(network: &LaminarNetwork<M>) -> Self {
        Self {
            node_flow: vec![0.0; network.topology.node_count],
            node_scale: vec![0.0; network.topology.node_count],
            component_flow: vec![0.0; network.topology.pressure_nodes.len()],
            component_scale: vec![0.0; network.topology.pressure_nodes.len()],
        }
    }

    pub(super) fn audit<M: FluidTransport>(
        &mut self,
        network: &LaminarNetwork<M>,
        records: &[LaminarNetworkBranchStep],
        output: &[f64],
        average: bool,
        coordinate_bound: f64,
    ) -> Result<(f64, f64), LaminarPipeError> {
        self.node_flow.fill(0.0);
        self.node_scale.fill(0.0);
        self.component_flow.fill(0.0);
        self.component_scale.fill(0.0);
        for (edge, report) in records.iter().enumerate() {
            let flow = if average {
                report.average_volume_flow_m3_s
            } else {
                report.final_volume_flow_m3_s
            };
            let scale = flow.abs() + network.flow_norms[edge] * coordinate_bound;
            let [from, to] = network.topology.branch_nodes[edge];
            self.node_flow[from] += flow;
            self.node_flow[to] -= flow;
            self.node_scale[from] += scale;
            self.node_scale[to] += scale;
        }
        let tolerance = ROUNDOFF * network.system.state_dim() as f64;
        let mut internal = 0.0_f64;
        for node in 0..self.node_flow.len() {
            let value = self.node_flow[node];
            let scale = self.node_scale[node];
            if !value.is_finite() || !scale.is_finite() {
                return Err(LaminarPipeError::Invalid(
                    "network continuity account is unrepresentable",
                ));
            }
            if let Some(port) = network.basis.boundary_slot[node] {
                let boundary_output = if average { output[port] } else { value };
                if (boundary_output - value).abs() > tolerance * scale {
                    return Err(LaminarPipeError::Invalid(
                        "native reservoir flow disagrees with original branches",
                    ));
                }
                let root = network.basis.pressure_reference[port];
                self.component_flow[root] += boundary_output;
                self.component_scale[root] += scale;
            } else {
                internal = internal.max(value.abs());
                if value.abs() > tolerance * scale {
                    return Err(LaminarPipeError::Invalid(
                        "network trajectory violates internal volume continuity",
                    ));
                }
            }
        }
        let mut boundary = 0.0_f64;
        for (&value, &scale) in self.component_flow.iter().zip(&self.component_scale) {
            boundary = boundary.max(value.abs());
            if !value.is_finite() || !scale.is_finite() || value.abs() > tolerance * scale {
                return Err(LaminarPipeError::Invalid(
                    "network component reservoir flows do not conserve volume",
                ));
            }
        }
        Ok((internal, boundary))
    }
}

pub(super) fn prepared_error(error: PreparedStepError) -> LaminarPipeError {
    match error {
        PreparedStepError::Solver(error) => LaminarPipeError::Native(error),
        PreparedStepError::Cancelled => {
            LaminarPipeError::Execution(super::super::GraphError::Interrupted)
        }
    }
}
