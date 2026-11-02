//! Fixed complete-network force records, using the original exact/heat owners.
use super::*;
use fs_phs::linear_exact::{
    DenseExponentialBudget, PreparedSourceLoss, prepared::PreparedQuadraticFlow,
};

pub(crate) struct PreparedLinearForcing {
    block: ExactBlock,
    flow: PreparedQuadraticFlow,
    component: Option<PreparedSourceLoss>,
    links: Vec<Option<PreparedSourceLoss>>,
    material: Option<PreparedSourceLoss>,
    radiation: Option<PreparedSourceLoss>,
    physical: usize,
    weights: Vec<f64>,
    rest_scale: f64,
    state: Vec<f64>,
    input: Vec<f64>,
    z: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: [f64; 1],
    heat_coordinates: Vec<f64>,
}

fn prepared_error(error: fs_phs::PreparedStepError) -> ModalCouplingError {
    match error {
        fs_phs::PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
        fs_phs::PreparedStepError::Solver(error) => ModalCouplingError::Propagation(error),
    }
}

impl PreparedLinearForcing {
    pub(crate) fn new(
        network: &CoupledModalSystem,
        budget: DenseExponentialBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        network.require_time_domain()?;
        let count = network.mode_count();
        let physical = count
            .checked_mul(2)
            .and_then(|v| v.checked_add(network.material_memory.len()))
            .and_then(|v| v.checked_add(network.radiation_state.len()))
            .ok_or_else(|| invalid("complete ramp physical dimension overflow"))?;
        // p_i = dt*F_i(t), c = dt*max|rest_force|. Impulse coordinates avoid
        // force derivatives proportional to 1/dt and keep the static preload
        // outside the external-work form, even for a nearly balanced preload.
        let n = physical
            .checked_add(count)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| invalid("complete force-record dimension overflow"))?;
        let d = n
            .checked_mul(2)
            .ok_or_else(|| invalid("complete positive ramp heat dimension overflow"))?;
        let channels = 3usize
            .checked_add(
                network
                    .connections
                    .iter()
                    .filter(|link| link.damping_n_s_m != 0.0)
                    .count(),
            )
            .ok_or_else(|| invalid("complete ramp source count overflow"))?;
        // Reserve every cache and original assembly, plus each native owner's
        // conservative workspace, BEFORE allocating any complete dense matrix.
        let bytes = d
            .checked_mul(d)
            .and_then(|v| v.checked_mul(8))
            .and_then(|v| v.checked_mul(channels.checked_mul(40)?.checked_add(64)?));
        let work = (d as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(4096))
            .and_then(|v| v.checked_mul(channels as u64 + 2));
        if d > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || d > budget.max_dimension
            || bytes.is_none_or(|v| v > budget.max_dense_bytes)
            || work.is_none_or(|v| v > budget.max_multiply_terms)
        {
            return Err(invalid(
                "complete force ramp and positive heat exceed explicit resources; no states were dropped",
            ));
        }
        let native_budget = DenseExponentialBudget {
            max_dimension: budget.max_dimension,
            max_dense_bytes: budget.max_dense_bytes / (channels + 2),
            max_multiply_terms: budget.max_multiply_terms / (channels as u64 + 2),
        };
        let modes: Vec<_> = network
            .models
            .iter()
            .flat_map(|m| m.modes().iter().copied())
            .collect();
        let damping: Vec<_> = (0..network.connections.len())
            .map(|i| network.connection_damping_column(i).map(<[f64]>::to_vec))
            .collect();
        let poles: Vec<_> = network.radiation.as_ref().map_or_else(Vec::new, |air| {
            air.model
                .poles
                .iter()
                .map(|pole| {
                    let mut coupling = vec![0.0; count];
                    for (&i, &value) in air.source_modes.iter().zip(&pole.coupling) {
                        coupling[i] = value;
                    }
                    crate::render::plate::impact::radiation::Pole {
                        omega: pole.omega,
                        zeta: pole.zeta,
                        coupling,
                    }
                })
                .collect()
        });
        let block = ExactBlock::new_with_damping_columns(
            &modes,
            &network.columns,
            &network.connections,
            network.relaxation_branches.clone(),
            &poles,
            Some(&damping),
            gate,
        )?;
        if block.system.state_dim() != physical {
            return Err(invalid(
                "force ramp lost complete material or acoustic state",
            ));
        }
        let (j, r, g) = block.system.structure();
        let q = block.split.quadratic_matrix();
        let inv_dt = finite(1.0 / network.dt)?;
        let rest_scale = block
            .rest_force
            .iter()
            .map(|f| f.abs())
            .fold(0.0_f64, f64::max);
        let mut generator = vec![0.0; n * n];
        let mut input_map = vec![0.0; n * count];
        let mut work_form = vec![0.0; n * n];
        for i in 0..physical {
            poll(Some(gate))?;
            for k in 0..physical {
                let jr = j[i * physical + k] - r[i * physical + k];
                for column in 0..physical {
                    generator[i * n + column] =
                        finite(generator[i * n + column] + jr * q[k * physical + column])?;
                }
            }
            for p in 0..count {
                generator[i * n + physical + p] = finite(g[i * count + p] * inv_dt)?;
                if rest_scale > 0.0 {
                    generator[i * n + n - 1] = finite(
                        generator[i * n + n - 1]
                            + g[i * count + p] * (block.rest_force[p] / rest_scale) * inv_dt,
                    )?;
                }
            }
        }
        for p in 0..count {
            input_map[(physical + p) * count + p] = inv_dt;
            work_form[(2 * p + 1) * n + physical + p] = 0.5 * inv_dt;
            work_form[(physical + p) * n + 2 * p + 1] = 0.5 * inv_dt;
        }
        let flow = PreparedQuadraticFlow::new_extended_input_controlled(
            &block.system,
            &block.split,
            &generator,
            &vec![1.0; count + 1],
            &input_map,
            network.dt,
            &[work_form],
            &[],
            native_budget,
            || gate.is_requested(),
        )
        .map_err(prepared_error)?;

        // Balance source integration using original storage coordinate scales;
        // no eigenvalue clipping or reconstructed resistance supplies a source.
        let mut weights = vec![1.0; n];
        for i in 0..physical {
            if q[i * physical + i] > 0.0 {
                weights[i] = finite(det::sqrt(q[i * physical + i]))?;
            }
        }
        for i in 0..n {
            poll(Some(gate))?;
            for k in 0..n {
                generator[i * n + k] = finite(generator[i * n + k] * weights[i] / weights[k])?;
            }
            for p in 0..count {
                input_map[i * count + p] = finite(input_map[i * count + p] * weights[i])?;
            }
        }
        let source = |rows: Vec<f64>| -> Result<Option<PreparedSourceLoss>, ModalCouplingError> {
            if rows.is_empty() {
                return Ok(None);
            }
            PreparedSourceLoss::new_controlled(
                n,
                count,
                &generator,
                &input_map,
                &rows,
                network.dt,
                native_budget,
                || gate.is_requested(),
            )
            .map(Some)
            .map_err(prepared_error)
        };
        let diagonal_rows = |resistance: &[f64]| -> Result<Vec<f64>, ModalCouplingError> {
            let mut rows = Vec::new();
            if resistance.is_empty() {
                return Ok(rows);
            }
            for i in 0..physical {
                let value = resistance[i * physical + i];
                if value == 0.0 {
                    continue;
                }
                let root = finite(det::sqrt(value))?;
                for k in 0..n {
                    rows.push(if k < physical {
                        finite(root * q[i * physical + k] / weights[k])?
                    } else {
                        0.0
                    });
                }
            }
            Ok(rows)
        };
        let component = source(diagonal_rows(&block.component_r)?)?;
        let material = source(diagonal_rows(&block.material_r)?)?;
        let radiation = source(diagonal_rows(&block.radiation_r)?)?;
        let mut links = Vec::with_capacity(block.link_loss_rows.len());
        for row in &block.link_loss_rows {
            poll(Some(gate))?;
            let mut pulled = Vec::new();
            if let Some(row) = row {
                pulled.resize(n, 0.0);
                for k in 0..physical {
                    pulled[k] = finite(
                        (0..physical)
                            .map(|i| row[i] * q[i * physical + k])
                            .sum::<f64>()
                            / weights[k],
                    )?;
                }
            }
            links.push(source(pulled)?);
        }
        Ok(Self {
            block,
            flow,
            component,
            links,
            material,
            radiation,
            physical,
            weights,
            rest_scale,
            state: vec![0.0; n],
            input: vec![0.0; count],
            z: vec![0.0; n + count + 1],
            balanced: vec![0.0; n],
            endpoint: vec![0.0; n],
            integrals: [0.0],
            heat_coordinates: vec![0.0; 2 * n],
        })
    }

    pub(crate) fn stage(
        &mut self,
        network: &mut CoupledModalSystem,
        start: &[f64],
        end: &[f64],
        gate: &CancelGate,
    ) -> Result<(), ModalCouplingError> {
        poll(Some(gate))?;
        let count = network.mode_count();
        if start.len() != count
            || end.len() != count
            || start.iter().chain(end).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "force ramp endpoints must cover every original finite modal force",
            ));
        }
        let x = flatten(&network.models);
        self.state[..2 * count].copy_from_slice(&x);
        let memory_end = 2 * count + network.material_memory.len();
        self.state[2 * count..memory_end].copy_from_slice(&network.material_memory);
        self.state[memory_end..self.physical].copy_from_slice(&network.radiation_state);
        let mut envelope = Vec::with_capacity(count);
        for p in 0..count {
            self.state[self.physical + p] = finite(network.dt * start[p])?;
            let delta = end[p] - start[p];
            self.input[p] = finite(if delta.is_finite() {
                network.dt * delta
            } else {
                network.dt * end[p] - network.dt * start[p]
            })?;
            envelope.push(
                finite(start[p] + self.block.rest_force[p])?
                    .abs()
                    .max(finite(end[p] + self.block.rest_force[p])?.abs()),
            );
        }
        *self.state.last_mut().expect("admitted full force state") =
            finite(network.dt * self.rest_scale)?;
        self.flow
            .evaluate_into(
                &self.state,
                &self.input,
                &mut self.z,
                &mut self.balanced,
                &mut self.endpoint,
                &mut self.integrals,
            )
            .map_err(ModalCouplingError::Propagation)?;
        for (out, (&value, &weight)) in self
            .balanced
            .iter_mut()
            .zip(self.state.iter().zip(&self.weights))
        {
            *out = finite(value * weight)?;
        }
        let mut heat = |channel: &Option<PreparedSourceLoss>| -> Result<f64, ModalCouplingError> {
            match channel {
                Some(owner) => owner
                    .evaluate_into_controlled(
                        &self.balanced,
                        &self.input,
                        &mut self.heat_coordinates,
                        || gate.is_requested(),
                    )
                    .map_err(prepared_error),
                None => Ok(0.0),
            }
        };
        let component = heat(&self.component)?;
        let material = heat(&self.material)?;
        let radiation = heat(&self.radiation)?;
        let mut connection = 0.0;
        let mut transfers = Vec::with_capacity(self.links.len());
        for (i, owner) in self.links.iter().enumerate() {
            poll(Some(gate))?;
            let loss = heat(owner)?;
            connection = finite(connection + loss)?;
            let link = &network.connections[i];
            let column = &network.columns[i];
            let before = finite(
                column
                    .iter()
                    .enumerate()
                    .map(|(k, b)| b * self.state[2 * k])
                    .sum::<f64>()
                    - link.rest_extension_m,
            )?;
            let after = finite(
                column
                    .iter()
                    .enumerate()
                    .map(|(k, b)| b * self.endpoint[2 * k])
                    .sum::<f64>()
                    - link.rest_extension_m,
            )?;
            let work =
                finite(-0.5 * link.stiffness_n_m * (after - before) * (after + before) - loss)?;
            transfers.push(crate::port_adoption::PortTransfer::work(work, false));
        }
        let candidate = PropagatedCandidate {
            endpoint: self.endpoint[..self.physical].to_vec(),
            start: self.state[..self.physical].to_vec(),
            input: envelope,
            losses: [component, connection, material, radiation],
            initial_bound_energy: finite(
                self.block.system.hamiltonian(&self.state[..self.physical]),
            )?,
            external_work: Some(self.integrals[0]),
            connection_transfers: Some(transfers),
        };
        // stage_candidate deliberately does not publish a held-input material
        // trajectory for this varying force. All physical history, source-domain
        // admission and work/heat checks still use this actual ramp endpoint.
        stage_candidate(network, end, candidate, Some(gate))
    }
}
