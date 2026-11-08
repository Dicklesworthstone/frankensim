//! Continuous reservoir and pump pressures on the original constrained network.
use super::super::tones::LaminarPressureTone;
use super::prepared::{NetworkFlowAudit, prepared_error};
use super::*;
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, PreparedExactTones};

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Prepare continuous cosine/sine pressure sources on every original port.
    /// Reservoir ports come first, followed by one directed pressure head per
    /// branch, exactly as in `system()`. Coefficients are port-major and every
    /// port supplies all prepared frequencies, including explicit zero heads.
    ///
    /// Every positive frequency must fit every original wall-memory band. All
    /// original profile states, positive heat channels and topology constraints
    /// remain. The explicit dense budget includes adapter scratch and channel
    /// preparation; no modes, pipes or forcing quadratures are discarded.
    pub fn prepare_tones(
        &self,
        cx: &Cx<'_>,
        frequencies_rad_s: &[f64],
        duration_s: f64,
        budget: DenseExponentialBudget,
    ) -> Result<PreparedLaminarNetworkTones<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        let (n, ports, edges) = (
            self.system.state_dim(),
            self.port_count(),
            self.branches.len(),
        );
        let coefficients =
            ports
                .checked_mul(frequencies_rad_s.len())
                .ok_or(LaminarPipeError::Invalid(
                    "network pressure-tone extent overflow",
                ))?;
        if !positive(duration_s)
            || frequencies_rad_s.is_empty()
            || coefficients > MAX_LAMINAR_NETWORK_AUGMENTED
            || frequencies_rad_s.iter().any(|w| !w.is_finite() || *w < 0.0)
            || frequencies_rad_s.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(LaminarPipeError::Invalid(
                "bounded complete network tones, increasing nonnegative frequencies and positive interval required",
            ));
        }
        // Reserve conservative channel matrices/factors and all transactional
        // buffers BEFORE constructing them. The native owner admits its own
        // complete exponential/power caches against the remaining budget.
        let extra_bytes = 256
            * (n + ports + edges + coefficients + self.topology.node_count) as u128
            + 64 * edges as u128 * n as u128 * n as u128;
        let extra_work = 128 * edges as u128 * (n as u128).pow(3)
            + 64 * n as u128 * (n + ports + coefficients) as u128;
        if extra_bytes > budget.max_dense_bytes as u128
            || extra_work > budget.max_multiply_terms as u128
        {
            return Err(LaminarPipeError::Invalid(
                "network tone channels exceed explicit capacity",
            ));
        }
        let native_budget = DenseExponentialBudget {
            max_dimension: budget.max_dimension,
            max_dense_bytes: budget.max_dense_bytes - extra_bytes as usize,
            max_multiply_terms: budget.max_multiply_terms - extra_work as u64,
        };
        for &omega in frequencies_rad_s {
            for pipe in &self.branches {
                checkpoint(cx)?;
                if omega > 0.0 {
                    pipe.impedance_at(omega)?;
                }
            }
        }
        let mut channels = Vec::with_capacity(edges);
        for source in &self.source_rows {
            checkpoint(cx)?;
            let mut matrix = vec![0.0; n * n];
            for row in source.chunks_exact(n) {
                checkpoint(cx)?;
                for i in 0..n {
                    for j in 0..=i {
                        let value = row[i] * row[j];
                        matrix[i * n + j] += value;
                        if i != j {
                            matrix[j * n + i] += value;
                        }
                    }
                }
            }
            channels.push(ResistanceWitness::try_from_source_rows(
                n,
                matrix,
                source.clone(),
            )?);
        }
        let native = PreparedExactTones::new_with_resistance_channels_controlled(
            &self.system,
            frequencies_rad_s,
            duration_s,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: 1.0,
            },
            native_budget,
            &channels,
            || cx.checkpoint().is_err(),
        )
        .map_err(prepared_error)?;
        let prepared = PreparedLaminarNetworkTones {
            network: self,
            native,
            cosine: vec![0.0; coefficients],
            sine: vec![0.0; coefficients],
            pending_pressure: vec![
                LaminarPressureTone {
                    cosine: Pressure::new(0.0),
                    sine: Pressure::new(0.0)
                };
                coefficients
            ],
            pressure_values: vec![Pressure::new(0.0); ports],
            input: vec![0.0; ports],
            dc: vec![0.0; ports],
            output: vec![0.0; ports],
            pending: vec![0.0; n],
            steady: vec![0.0; n],
            drive: vec![0.0; n],
            difference: vec![0.0; n],
            residual: vec![0.0; n],
            heat: vec![0.0; edges],
            pending_work: vec![0.0; ports],
            work: vec![0.0; ports],
            records: vec![LaminarNetworkBranchStep::default(); edges],
            pending_records: vec![LaminarNetworkBranchStep::default(); edges],
            boundary_flows: vec![0.0; self.topology.pressure_nodes.len()],
            flow_audit: NetworkFlowAudit::new(self),
            last: None,
        };
        checkpoint(cx)?;
        Ok(prepared)
    }
}

/// One complete coupled trajectory with per-pipe heat and per-source work.
/// Preparation owns all scratch; successful warm calls reuse every buffer.
pub struct PreparedLaminarNetworkTones<'a, M> {
    network: &'a LaminarNetwork<M>,
    native: PreparedExactTones<'a>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    pending_pressure: Vec<LaminarPressureTone>,
    pressure_values: Vec<Pressure>,
    input: Vec<f64>,
    dc: Vec<f64>,
    output: Vec<f64>,
    pending: Vec<f64>,
    steady: Vec<f64>,
    drive: Vec<f64>,
    difference: Vec<f64>,
    residual: Vec<f64>,
    heat: Vec<f64>,
    pending_work: Vec<f64>,
    work: Vec<f64>,
    records: Vec<LaminarNetworkBranchStep>,
    pending_records: Vec<LaminarNetworkBranchStep>,
    boundary_flows: Vec<f64>,
    flow_audit: NetworkFlowAudit,
    last: Option<LaminarNetworkStep>,
}

impl<M: FluidTransport> PreparedLaminarNetworkTones<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.native.angular_frequencies_rad_s()
    }
    pub fn port_count(&self) -> usize {
        self.native.port_dim()
    }
    pub fn branch_steps(&self) -> &[LaminarNetworkBranchStep] {
        if self.last.is_some() {
            &self.records
        } else {
            &[]
        }
    }
    pub fn boundary_volume_flows_m3_s(&self) -> &[f64] {
        if self.last.is_some() {
            &self.boundary_flows
        } else {
            &[]
        }
    }
    /// Signed work in original port order. Reservoir pressures use each
    /// component's original reference gauge; directed head work follows it.
    /// The total reservoir work is gauge-invariant, individual ports need not be.
    pub fn port_work_j(&self) -> &[f64] {
        if self.last.is_some() { &self.work } else { &[] }
    }
    pub const fn last_step(&self) -> Option<LaminarNetworkStep> {
        self.last
    }

    /// Advance the actual continuous pressures and publish their end phases
    /// together with accepted state, clock, branch records and source accounts.
    /// The flat source buffer has `port_count()*tone_count` entries. Each
    /// interval starts at the accepted state's clock; retain the rotated buffer
    /// for continuous continuation. A DC tone requires exactly zero sine.
    ///
    /// Passive contraction bounds the entire trajectory, including pressure
    /// peaks between endpoints. The DC equilibrium improves this bound for
    /// steady bias. Oscillatory amplitudes use a conservative forcing-norm
    /// envelope, which can refuse long intervals or strongly cancelling tones.
    /// This is exact forcing integration within the original finite wall model,
    /// not a new radial-continuum or turbulence certificate.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarNetworkState,
        pressures: &mut [LaminarPressureTone],
    ) -> Result<LaminarNetworkStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.network.admit_owner(state)?;
        if pressures.len() != self.cosine.len() {
            return Err(LaminarPipeError::Invalid(
                "complete port-major network pressure tones required",
            ));
        }
        let dt = self.duration_s();
        let next_time = state.accepted_time_s + dt;
        if !next_time.is_finite() || next_time <= state.accepted_time_s {
            return Err(LaminarPipeError::Invalid(
                "network physical clock must advance finitely",
            ));
        }
        self.pack(cx, pressures)?;
        let coordinate_bound = self.admit_interval(cx, state)?;
        let native = self
            .native
            .step_accounted_into_controlled(
                &state.coordinates,
                &self.cosine,
                &self.sine,
                &mut self.pending,
                &mut self.output,
                &mut self.heat,
                &mut self.pending_work,
                || cx.checkpoint().is_err(),
            )
            .map_err(prepared_error)?;
        let n = self.pending.len();
        let boundary = self.network.topology.pressure_nodes.len();
        for edge in 0..self.records.len() {
            checkpoint(cx)?;
            let row = &self.network.flow_rows[edge];
            let r = self.network.basis.columns;
            let mean = &self.network.basis.energy_rows[edge * r..(edge + 1) * r];
            let initial_mean = dot(mean, &state.coordinates[..r]);
            let final_mean = dot(mean, &self.pending[..r]);
            let start = self.network.offsets[edge];
            let end = self.network.offsets[edge + 1];
            let record = &mut self.pending_records[edge];
            record.initial_volume_flow_m3_s = dot(row, &state.coordinates);
            record.final_volume_flow_m3_s = dot(row, &self.pending);
            record.average_volume_flow_m3_s = self.output[boundary + edge];
            record.displaced_volume_m3 = dt * record.average_volume_flow_m3_s;
            record.initial_mean_kinetic_energy_j = 0.5 * initial_mean * initial_mean;
            record.final_mean_kinetic_energy_j = 0.5 * final_mean * final_mean;
            record.initial_profile_memory_energy_j = energy(&state.coordinates[start..end]);
            record.final_profile_memory_energy_j = energy(&self.pending[start..end]);
            record.viscous_heat_j = self.heat[edge];
            if ![
                record.initial_volume_flow_m3_s,
                record.final_volume_flow_m3_s,
                record.average_volume_flow_m3_s,
                record.displaced_volume_m3,
                record.initial_mean_kinetic_energy_j,
                record.final_mean_kinetic_energy_j,
                record.initial_profile_memory_energy_j,
                record.final_profile_memory_energy_j,
                record.viscous_heat_j,
            ]
            .into_iter()
            .all(f64::is_finite)
                || record.viscous_heat_j < 0.0
                || (record.average_volume_flow_m3_s != 0.0 && record.displaced_volume_m3 == 0.0)
            {
                return Err(LaminarPipeError::Invalid(
                    "network tone branch flow, heat or storage is unrepresentable",
                ));
            }
        }
        let local_heat = accurate_sum(self.heat.iter().copied());
        if !local_heat.is_finite()
            || (local_heat - native.dissipated).abs()
                > ROUNDOFF * n as f64 * local_heat.abs().max(native.dissipated.abs())
        {
            return Err(LaminarPipeError::Invalid(
                "original branch heat does not close continuous network account",
            ));
        }
        let endpoint = self.flow_audit.audit(
            self.network,
            &self.pending_records,
            &self.output,
            false,
            coordinate_bound,
        )?;
        let average = self.flow_audit.audit(
            self.network,
            &self.pending_records,
            &self.output,
            true,
            coordinate_bound,
        )?;
        let reservoir_work_j = accurate_sum(self.pending_work[..boundary].iter().copied());
        let branch_head_work_j = accurate_sum(self.pending_work[boundary..].iter().copied());
        let work_scale =
            accurate_sum(self.pending_work.iter().map(|v| v.abs())) + native.supplied.abs();
        if ![reservoir_work_j, branch_head_work_j, work_scale]
            .into_iter()
            .all(f64::is_finite)
            || (reservoir_work_j + branch_head_work_j - native.supplied).abs()
                > ROUNDOFF * n as f64 * work_scale
        {
            return Err(LaminarPipeError::Invalid(
                "continuous reservoir and head work do not close native supply",
            ));
        }
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
        checkpoint(cx)?;
        state.coordinates.copy_from_slice(&self.pending);
        state.accepted_time_s = next_time;
        pressures.copy_from_slice(&self.pending_pressure);
        core::mem::swap(&mut self.records, &mut self.pending_records);
        core::mem::swap(&mut self.work, &mut self.pending_work);
        self.boundary_flows
            .copy_from_slice(&self.output[..boundary]);
        self.last = Some(report);
        Ok(report)
    }

    fn pack(
        &mut self,
        cx: &Cx<'_>,
        pressures: &[LaminarPressureTone],
    ) -> Result<(), LaminarPipeError> {
        let count = self.native.tone_count();
        let boundary = self.network.topology.pressure_nodes.len();
        for tone in 0..count {
            checkpoint(cx)?;
            for port in 0..self.input.len() {
                let index = port * count + tone;
                let (c, s) = (
                    pressures[index].cosine.value(),
                    pressures[index].sine.value(),
                );
                if !c.is_finite()
                    || !s.is_finite()
                    || (self.native.angular_frequencies_rad_s()[tone] == 0.0 && s != 0.0)
                {
                    return Err(LaminarPipeError::Invalid(
                        "finite network tones and zero DC sine required",
                    ));
                }
                let (next_c, next_s) = self.native.advance_force_phase(tone, c, s)?;
                self.pending_pressure[index] = LaminarPressureTone {
                    cosine: Pressure::new(next_c),
                    sine: Pressure::new(next_s),
                };
                self.pressure_values[port] = pressures[index].cosine;
            }
            self.network.fill_inputs(
                &self.pressure_values[..boundary],
                Some(&self.pressure_values[boundary..]),
                &mut self.input,
            )?;
            for port in 0..self.input.len() {
                self.cosine[port * count + tone] = self.input[port];
            }
            for port in 0..self.input.len() {
                self.pressure_values[port] = pressures[port * count + tone].sine;
            }
            self.network.fill_inputs(
                &self.pressure_values[..boundary],
                Some(&self.pressure_values[boundary..]),
                &mut self.input,
            )?;
            for port in 0..self.input.len() {
                self.sine[port * count + tone] = self.input[port];
            }
        }
        Ok(())
    }

    fn admit_interval(
        &mut self,
        cx: &Cx<'_>,
        state: &LaminarNetworkState,
    ) -> Result<f64, LaminarPipeError> {
        let n = self.pending.len();
        let ports = self.input.len();
        let count = self.native.tone_count();
        let (_, resistance, g) = self.network.system.structure();
        self.dc.fill(0.0);
        let mut oscillatory = 0.0;
        for tone in 0..count {
            checkpoint(cx)?;
            if self.native.angular_frequencies_rad_s()[tone] == 0.0 {
                for port in 0..ports {
                    self.dc[port] = self.cosine[port * count + tone];
                }
                continue;
            }
            for i in 0..n {
                let row = &g[i * ports..(i + 1) * ports];
                let c = accurate_sum((0..ports).map(|p| row[p] * self.cosine[p * count + tone]));
                let s = accurate_sum((0..ports).map(|p| row[p] * self.sine[p * count + tone]));
                let cm = accurate_sum(
                    (0..ports).map(|p| (row[p] * self.cosine[p * count + tone]).abs()),
                );
                let sm =
                    accurate_sum((0..ports).map(|p| (row[p] * self.sine[p * count + tone]).abs()));
                self.residual[i] = det::hypot(
                    c.abs() + ROUNDOFF * ports as f64 * cm,
                    s.abs() + ROUNDOFF * ports as f64 * sm,
                );
            }
            oscillatory += norm(&self.residual);
        }
        for i in 0..n {
            self.steady[i] = dot(
                &self.network.steady_response[i * ports..(i + 1) * ports],
                &self.dc,
            );
            let row = &g[i * ports..(i + 1) * ports];
            self.drive[i] = dot(row, &self.dc);
            let magnitude = accurate_sum(row.iter().zip(&self.dc).map(|(g, u)| (g * u).abs()));
            self.residual[i] = self.drive[i].abs() + ROUNDOFF * ports as f64 * magnitude;
            self.difference[i] = state.coordinates[i] - self.steady[i];
        }
        let dt = self.duration_s();
        let transient = norm(&state.coordinates) + dt * (norm(&self.residual) + oscillatory);
        for i in 0..n {
            let row = &resistance[i * n..(i + 1) * n];
            let restoring = dot(row, &self.steady);
            let magnitude = accurate_sum(row.iter().zip(&self.steady).map(|(r, x)| (r * x).abs()));
            let drive_magnitude = accurate_sum(
                g[i * ports..(i + 1) * ports]
                    .iter()
                    .zip(&self.dc)
                    .map(|(g, u)| (g * u).abs()),
            );
            self.residual[i] = (self.drive[i] - restoring).abs()
                + ROUNDOFF * (n as f64 * magnitude + ports as f64 * drive_magnitude);
        }
        let deviation = norm(&self.difference) + dt * (norm(&self.residual) + oscillatory);
        let inflation = 1.0 + ROUNDOFF * (n + count) as f64;
        let coordinate_bound = transient.min(norm(&self.steady) + deviation) * inflation;
        if ![oscillatory, transient, deviation, coordinate_bound]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(LaminarPipeError::Invalid(
                "continuous network momentum envelope is unrepresentable",
            ));
        }
        for edge in 0..self.records.len() {
            checkpoint(cx)?;
            let row = &self.network.flow_rows[edge];
            let magnitude = accurate_sum(row.iter().zip(&self.steady).map(|(q, x)| (q * x).abs()));
            let observed = dot(row, &self.steady).abs() + ROUNDOFF * n as f64 * magnitude;
            let row_norm = self.network.flow_norms[edge];
            let bound = (row_norm * transient).min(observed + row_norm * deviation) * inflation;
            if !bound.is_finite() {
                return Err(LaminarPipeError::Invalid(
                    "network tone Reynolds envelope is unrepresentable",
                ));
            }
            self.pending_records[edge].reynolds_envelope =
                self.network.branches[edge].admit_flow(bound)?;
        }
        Ok(coordinate_bound)
    }
}

fn accurate_sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
mod tests;
