//! Full retained bilateral pencil; nonlinear contacts keep their existing
//! held discrete-gradient reactions around this SAME affine response.
use super::*;
use fs_math::det;
use fs_phs::linear_exact::{
    LinearSplit, LinearSplitBudget, integrated_resistance_loss, split_quadratic, step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepRecord};

pub(super) mod prepared;

pub(super) struct ExactNetwork {
    blocks: Vec<(Vec<usize>, ExactBlock)>,
    prepared: Option<prepared::PreparedNetwork>,
    count: usize,
}

// Diagonal component dynamics and held inputs cannot connect distinct K/C
// islands. A permutation therefore makes the full generator block diagonal;
// its exponential and each resistance integral are the direct sum of these
// exact blocks. Contacts still act through the complete original force vector.
// This API's original modes have identity mass in normalized coordinates;
// connection storage is k*(column.q-rest)^2/2 with no momentum cross terms.
// Thus these actual K/C columns exhaust every possible inter-mode coupling.
// This changes no retained coordinate, physical limit or dense-block capacity.
impl ExactNetwork {
    pub(super) fn new_prepared(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        count: usize,
        dt: f64,
        budget: PreparedModalBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let prepared = prepared::PreparedNetwork::new(
            models, columns, links, count, dt, budget, gate,
        )?;
        Ok(Self { blocks: Vec::new(), prepared: Some(prepared), count })
    }

    pub(super) fn new(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        count: usize,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let modes: Vec<_> = models
            .iter()
            .flat_map(|m| m.modes().iter().copied())
            .collect();
        if count == 0
            || modes.len() != count
            || columns.len() != links.len()
            || columns.iter().any(|c| c.len() != count)
        {
            return Err(invalid("complete modal island layout"));
        }
        let mut labels: Vec<_> = (0..count).collect();
        let root = |parents: &[usize], mut i: usize| {
            while parents[i] != i {
                i = parents[i];
            }
            i
        };
        for (column, link) in columns.iter().zip(links) {
            poll(Some(gate))?;
            if link.stiffness_n_m == 0.0 && link.damping_n_s_m == 0.0 {
                continue;
            }
            let mut members = column.iter().enumerate().filter(|(_, v)| **v != 0.0);
            if let Some((first, _)) = members.next() {
                for (next, _) in members {
                    let (a, b) = (root(&labels, first), root(&labels, next));
                    if a != b {
                        labels[a.max(b)] = a.min(b);
                    }
                }
            }
        }
        let mut groups = std::collections::BTreeMap::<usize, Vec<usize>>::new();
        for index in 0..count {
            groups.entry(root(&labels, index)).or_default().push(index);
        }
        let mut blocks = Vec::with_capacity(groups.len());
        for indices in groups.into_values() {
            poll(Some(gate))?;
            let local_modes: Vec<_> = indices.iter().map(|&i| modes[i]).collect();
            let local_columns: Vec<Vec<_>> = columns
                .iter()
                .map(|c| indices.iter().map(|&i| c[i]).collect())
                .collect();
            let block = ExactBlock::new(&local_modes, &local_columns, links, gate)?;
            blocks.push((indices, block));
        }
        Ok(Self { blocks, prepared: None, count })
    }

    fn local_state(indices: &[usize], x: &[f64]) -> Vec<f64> {
        indices
            .iter()
            .flat_map(|&i| [x[2 * i], x[2 * i + 1]])
            .collect()
    }

    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if let Some(prepared) = &self.prepared {
            return prepared.hamiltonian(x);
        }
        self.blocks
            .iter()
            .map(|(indices, b)| b.system.hamiltonian(&Self::local_state(indices, x)))
            .sum()
    }

    pub(super) fn advance(
        &self,
        models: &[ModalAcousticTimeModel],
        external: &[f64],
        dt: f64,
    ) -> Result<(StepRecord, Vec<f64>, Vec<f64>), ModalCouplingError> {
        if external.len() != self.count || external.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "external forces must match every original modal coordinate",
            ));
        }
        if let Some(prepared) = &self.prepared {
            return prepared.advance(models, external, dt);
        }
        let x = flatten(models);
        let mut u = vec![0.0; self.count];
        let mut record = StepRecord {
            x: vec![0.0; 2 * self.count],
            y: vec![0.0; self.count],
            delta_h: 0.0,
            dissipated: 0.0,
            supplied: 0.0,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        };
        for (indices, block) in &self.blocks {
            let local_external: Vec<_> = indices.iter().map(|&i| external[i]).collect();
            let local_u = block.held_input(&local_external)?;
            let local = step_render(
                &block.system,
                &block.split,
                &Self::local_state(indices, &x),
                &local_u,
                dt,
            )
            .map_err(ModalCouplingError::Propagation)?;
            for (j, &i) in indices.iter().enumerate() {
                record.x[2 * i..2 * i + 2].copy_from_slice(&local.x[2 * j..2 * j + 2]);
                record.y[i] = local.y[j];
                u[i] = local_u[j];
            }
            record.delta_h = finite(record.delta_h + local.delta_h)?;
            record.dissipated = finite(record.dissipated + local.dissipated)?;
            record.supplied = finite(record.supplied + local.supplied)?;
        }
        Ok((record, x, u))
    }

    /// The free prediction of the next sample: (end state, start state). With the prepared owner and `verify == false` the
    /// loss Gramians and the trajectory energy-residual check are skipped (see `PreparedNetwork::predict_free`); `verify`
    /// runs the full `advance`. Every other owner always runs the full `advance`.
    pub(super) fn predict(&self, models: &[ModalAcousticTimeModel], external: &[f64], dt: f64, verify: bool)
        -> Result<(Vec<f64>, Vec<f64>), ModalCouplingError>
    {
        if let (Some(prepared), false) = (&self.prepared, verify) {
            if external.len() != self.count || external.iter().any(|v| !v.is_finite()) {
                return Err(invalid("external forces must match every original modal coordinate"));
            }
            return prepared.predict_free(models, external, dt);
        }
        let (record, x, _) = self.advance(models, external, dt)?;
        Ok((record.x, x))
    }

    pub(super) fn response(&self, column: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        if column.len() != self.count {
            return Err(invalid("complete contact force column"));
        }
        if let Some(prepared) = &self.prepared {
            return prepared.response(column, dt);
        }
        let mut response = vec![0.0; self.count];
        for (indices, block) in &self.blocks {
            let local_column: Vec<_> = indices.iter().map(|&i| column[i]).collect();
            for (&i, value) in indices.iter().zip(block.response(&local_column, dt)?) {
                response[i] = value;
            }
        }
        Ok(response)
    }

    fn losses(
        &self,
        x: &[f64],
        u: &[f64],
        dt: f64,
        total: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        if let Some(prepared) = &self.prepared {
            return prepared.losses(x, u, dt, total);
        }
        let (mut component, mut connection) = (0.0, 0.0);
        for (indices, block) in &self.blocks {
            let local_x = Self::local_state(indices, x);
            let local_u: Vec<_> = indices.iter().map(|&i| u[i]).collect();
            let block_total = if self.blocks.len() == 1 {
                total
            } else {
                // Recover this island's independently integrated total for
                // the existing channel split, never apportion a global loss.
                step_render(&block.system, &block.split, &local_x, &local_u, dt)
                    .map_err(ModalCouplingError::Propagation)?
                    .dissipated
            };
            let (a, b) = block.losses(&local_x, &local_u, dt, block_total)?;
            component = finite(component + a)?;
            connection = finite(connection + b)?;
        }
        Ok((component, connection))
    }

    fn connection_work(
        &self,
        x: &[f64],
        u: &[f64],
        endpoint: &[f64],
        dt: f64,
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        connection_loss: f64,
        gate: Option<&CancelGate>,
    ) -> Result<Vec<crate::port_adoption::PortTransfer>, ModalCouplingError> {
        if let Some(prepared) = &self.prepared {
            return prepared.connection_work(x, u, endpoint, dt, columns, links, connection_loss, gate);
        }
        let mut work = vec![0.0; links.len()];
        for (indices, block) in &self.blocks {
            poll(gate)?;
            let local_x = Self::local_state(indices, x);
            let local_endpoint = Self::local_state(indices, endpoint);
            let local_u: Vec<_> = indices.iter().map(|&i| u[i]).collect();
            let local_columns: Vec<Vec<_>> = columns
                .iter()
                .map(|c| indices.iter().map(|&i| c[i]).collect())
                .collect();
            let local_loss = if self.blocks.len() == 1 {
                connection_loss
            } else if block.connection_r.iter().all(|v| *v == 0.0) {
                0.0
            } else {
                integrated_resistance_loss(
                    &block.system,
                    &block.split,
                    &local_x,
                    &local_u,
                    dt,
                    &block.connection_r,
                )
                .map_err(ModalCouplingError::Propagation)?
            };
            let transfers = block.connection_work(
                &local_x,
                &local_u,
                &local_endpoint,
                dt,
                &local_columns,
                links,
                local_loss,
                gate,
            )?;
            for (sum, transfer) in work.iter_mut().zip(transfers) {
                *sum = finite(*sum + transfer.energy_j)?;
            }
        }
        Ok(work
            .into_iter()
            .map(|value| crate::port_adoption::PortTransfer::work(value, false))
            .collect())
    }
}

struct ExactBlock {
    system: PortHamiltonian,
    split: LinearSplit,
    rest_force: Vec<f64>,
    component_r: Vec<f64>,
    connection_r: Vec<f64>,
}

impl ExactBlock {
    fn new(
        modes: &[crate::modal_acoustic_time::ModalAcousticMode],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let count = modes.len();
        let n = count
            .checked_mul(2)
            .ok_or_else(|| invalid("full linear state extent"))?;
        let block = n
            .checked_add(count)
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(2))
            .ok_or_else(|| invalid("full linear exponential block overflow"))?;
        if block > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(ModalCouplingError::Propagation(
                fs_phs::PhsError::LinearSplitBudget {
                    what: "complete modal network exceeds dense exponential block capacity; no modes were dropped",
                },
            ));
        }
        let mut q = vec![0.0; n * n];
        let mut j = vec![0.0; n * n];
        let mut r = vec![0.0; n * n];
        let mut g = vec![0.0; n * count];
        let mut at = 0;
        for mode in modes {
            poll(Some(gate))?;
            let qi = 2 * at;
            let pi = qi + 1;
            q[qi * n + qi] = mode.angular_frequency_rad_s.powi(2);
            q[pi * n + pi] = 1.0;
            j[qi * n + pi] = 1.0;
            j[pi * n + qi] = -1.0;
            g[pi * count + at] = 1.0;
            r[pi * n + pi] = 2.0 * mode.damping_ratio * mode.angular_frequency_rad_s;
            at += 1;
        }
        let component_r = r.clone();
        let mut connection_r = vec![0.0; n * n];
        let mut rest_force = vec![0.0; count];
        for (column, link) in columns.iter().zip(links) {
            poll(Some(gate))?;
            for a in 0..count {
                rest_force[a] += link.stiffness_n_m * link.rest_extension_m * column[a];
                for b in 0..count {
                    q[(2 * a) * n + 2 * b] += link.stiffness_n_m * column[a] * column[b];
                    connection_r[(2 * a + 1) * n + 2 * b + 1] +=
                        link.damping_n_s_m * column[a] * column[b];
                }
            }
        }
        for (value, extra) in r.iter_mut().zip(&connection_r) {
            *value += extra;
        }
        let storage = QuadraticStorage::new(q, n).map_err(ModalCouplingError::Propagation)?;
        let system = PortHamiltonian::new(n, count, j, r, g, Box::new(storage))
            .map_err(ModalCouplingError::Propagation)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: core::f64::consts::PI,
            },
        )
        .map_err(ModalCouplingError::Propagation)?;
        poll(Some(gate))?;
        Ok(Self {
            system,
            split,
            rest_force,
            component_r,
            connection_r,
        })
    }

    fn held_input(&self, external: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        if external.len() != self.rest_force.len() || external.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "external forces must match the complete finite linear pencil",
            ));
        }
        external
            .iter()
            .zip(&self.rest_force)
            .map(|(f, r)| finite(f + r))
            .collect()
    }

    pub(super) fn response(&self, column: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        // A derivative with respect to held force has zero initial state and
        // zero rest load. No finite-difference material or trajectory probes.
        let rec = step_render(
            &self.system,
            &self.split,
            &vec![0.0; self.system.state_dim()],
            column,
            dt,
        )
        .map_err(ModalCouplingError::Propagation)?;
        Ok(rec.x.chunks_exact(2).map(|s| s[0]).collect())
    }

    fn losses(
        &self,
        x: &[f64],
        u: &[f64],
        dt: f64,
        total: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        if self.component_r.iter().all(|v| *v == 0.0) {
            return Ok((0.0, total));
        }
        if self.connection_r.iter().all(|v| *v == 0.0) {
            return Ok((total, 0.0));
        }
        let component =
            integrated_resistance_loss(&self.system, &self.split, x, u, dt, &self.component_r)
                .map_err(ModalCouplingError::Propagation)?;
        let connection =
            integrated_resistance_loss(&self.system, &self.split, x, u, dt, &self.connection_r)
                .map_err(ModalCouplingError::Propagation)?;
        let scale = self
            .system
            .hamiltonian(x)
            .abs()
            .max(total.abs())
            .max(f64::MIN_POSITIVE);
        if (component + connection - total).abs() > 1e-12 * scale {
            return Err(invalid(
                "independent linear loss channels do not sum to the complete loss",
            ));
        }
        Ok((component, connection))
    }

    // Work INTO the retained bodies from each signed connection column:
    // integral[(-k*e-c*e_dot)*e_dot] = -delta(k*e^2/2)-integral[c*e_dot^2].
    // Endpoint reaction times displacement is not this exact-trajectory work.
    fn connection_work(
        &self,
        x: &[f64],
        u: &[f64],
        endpoint: &[f64],
        dt: f64,
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        connection_loss: f64,
        gate: Option<&CancelGate>,
    ) -> Result<Vec<crate::port_adoption::PortTransfer>, ModalCouplingError> {
        let damped = links
            .iter()
            .filter(|link| link.damping_n_s_m != 0.0)
            .count();
        let n = self.system.state_dim();
        columns
            .iter()
            .zip(links)
            .map(|(column, link)| {
                poll(gate)?;
                let coordinate = |state: &[f64]| -> Result<f64, ModalCouplingError> {
                    let sum = column
                        .iter()
                        .enumerate()
                        .try_fold(0.0, |s, (i, b)| finite(s + b * state[2 * i]))?;
                    finite(sum - link.rest_extension_m)
                };
                let before = coordinate(x)?;
                let after = coordinate(endpoint)?;
                let loss = if link.damping_n_s_m == 0.0 {
                    0.0
                } else if damped == 1 {
                    // The existing independently integrated connection channel is
                    // already exactly this one damper; do not integrate it twice.
                    connection_loss
                } else {
                    let mut resistance = vec![0.0; n * n];
                    for (a, left) in column.iter().enumerate() {
                        poll(gate)?;
                        for (b, right) in column.iter().enumerate() {
                            resistance[(2 * a + 1) * n + 2 * b + 1] =
                                link.damping_n_s_m * left * right;
                        }
                    }
                    integrated_resistance_loss(&self.system, &self.split, x, u, dt, &resistance)
                        .map_err(ModalCouplingError::Propagation)?
                };
                let work =
                    finite(-0.5 * link.stiffness_n_m * (after - before) * (after + before) - loss)?;
                Ok(crate::port_adoption::PortTransfer::work(work, false))
            })
            .collect()
    }
}

fn flatten(models: &[ModalAcousticTimeModel]) -> Vec<f64> {
    models
        .iter()
        .flat_map(|m| {
            m.states()
                .iter()
                .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
        })
        .collect()
}

pub(super) fn stage(
    network: &mut CoupledModalSystem,
    external: &[f64],
    gate: Option<&CancelGate>,
) -> Result<(), ModalCouplingError> {
    poll(gate)?;
    let sample = network
        .frame
        .sample
        .checked_add(1)
        .ok_or_else(|| invalid("coupled sample clock overflow"))?;
    let before = network.total_energy_j()?;
    let (rec_x, x, u, component_loss, connection_loss, h_x) = if let Some(exact) = &network.exact {
        let (rec, x, u) = exact.advance(&network.models, external, network.dt)?;
        poll(gate)?;
        let (comp_loss, conn_loss) = exact.losses(&x, &u, network.dt, rec.dissipated)?;
        let h_x = exact.hamiltonian(&x);
        (rec.x, x, u, comp_loss, conn_loss, h_x)
    } else {
        if external.len() != network.mode_count() || external.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "external forces must match the complete finite linear pencil",
            ));
        }
        let x = flatten(&network.models);
        let u = external.to_vec();
        let mut rec_x = Vec::with_capacity(x.len());
        let mut comp_loss = 0.0_f64;
        let mut k = 0usize;
        for model in &network.models {
            poll(gate)?;
            for (mode_idx, (&mode, &state)) in model.modes().iter().zip(model.states()).enumerate()
            {
                let next = model.advance_exact_zoh_step(mode_idx, state, external[k]);
                rec_x.push(next.displacement_m_sqrt_kg);
                rec_x.push(next.velocity_m_sqrt_kg_per_s);
                let w = mode.angular_frequency_rad_s;
                let e_before = 0.5
                    * (state.velocity_m_sqrt_kg_per_s.powi(2)
                        + (w * state.displacement_m_sqrt_kg).powi(2));
                let e_after = 0.5
                    * (next.velocity_m_sqrt_kg_per_s.powi(2)
                        + (w * next.displacement_m_sqrt_kg).powi(2));
                let dq = next.displacement_m_sqrt_kg - state.displacement_m_sqrt_kg;
                comp_loss = finite(comp_loss + (e_before + external[k] * dq - e_after))?;
                k += 1;
            }
        }
        (rec_x, x, u, comp_loss, 0.0, before)
    };
    let mut pressure = 0.0;
    let mut modal_after = 0.0;
    let mut work = 0.0;
    for i in 0..network.models.len() {
        poll(gate)?;
        let range = network.offsets[i]..network.offsets[i + 1];
        let states: Vec<_> = range
            .clone()
            .map(|k| ModalAcousticState {
                displacement_m_sqrt_kg: rec_x[2 * k],
                velocity_m_sqrt_kg_per_s: rec_x[2 * k + 1],
            })
            .collect();
        network.candidates[i]
            .restore_states(&states)
            .map_err(|source| ModalCouplingError::Component {
                component: i,
                source,
            })?;
        let mut component_pressure = 0.0;
        for ((k, mode), state) in range.zip(network.models[i].modes()).zip(&states) {
            network.old_q[k] = x[2 * k];
            network.free_delta[k] = finite(rec_x[2 * k] - x[2 * k])?;
            work = finite(work + external[k] * network.free_delta[k])?;
            modal_after = finite(
                modal_after
                    + 0.5
                        * (state.velocity_m_sqrt_kg_per_s.powi(2)
                            + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg)
                                .powi(2)),
            )?;
            component_pressure = finite(
                component_pressure
                    + mode.pressure_per_modal_velocity.re * state.velocity_m_sqrt_kg_per_s
                    + mode.pressure_per_modal_velocity.im
                        * mode.angular_frequency_rad_s
                        * state.displacement_m_sqrt_kg,
            )?;
        }
        limit(
            "component observer pressure",
            component_pressure.abs(),
            network.models[i].budget().maximum_abs_pressure_pa,
        )?;
        pressure = finite(pressure + component_pressure)?;
    }
    let mut spring_after = 0.0;
    network.forces.copy_from_slice(external);
    // H_quadratic(t) <= (sqrt(H0)+dt*||u||/sqrt(2))^2 for PSD loss.
    // Each individual spring obeys k*(b.q)^2 <= 2*H_quadratic.
    // This bounds reaction magnitude throughout the step, not just endpoints.
    let root_twice_energy =
        finite(det::sqrt(2.0 * h_x) + network.dt * det::sqrt(u.iter().map(|f| f * f).sum()))?;
    for (index, (column, link)) in network.columns.iter().zip(&network.connections).enumerate() {
        poll(gate)?;
        let extension = extension(&network.candidates, column, link.rest_extension_m)?;
        spring_after = finite(spring_after + 0.5 * link.stiffness_n_m * extension * extension)?;
        let velocity: f64 = column
            .iter()
            .enumerate()
            .map(|(i, b)| b * rec_x[2 * i + 1])
            .sum();
        network.reactions[index] =
            finite(-link.stiffness_n_m * extension - link.damping_n_s_m * velocity)?;
        for (force, b) in network.forces.iter_mut().zip(column) {
            *force = finite(*force + b * network.reactions[index])?;
        }
        let force_bound = finite(
            (det::sqrt(link.stiffness_n_m)
                + link.damping_n_s_m * det::sqrt(column.iter().map(|b| b * b).sum()))
                * root_twice_energy
                + link.stiffness_n_m * link.rest_extension_m.abs(),
        )?;
        limit(
            "connection force over complete exact interval",
            force_bound,
            network.config.maximum_abs_connection_force_n,
        )?;
    }
    let after = finite(modal_after + spring_after)?;
    limit("total energy", after, network.config.maximum_total_energy_j)?;
    limit(
        "observer pressure",
        pressure.abs(),
        network.config.maximum_abs_pressure_pa,
    )?;
    let residual = finite(after - before + component_loss + connection_loss - work)?;
    let scale = before
        .max(after)
        .max(work.abs())
        .max((component_loss + connection_loss).abs());
    let tolerance = finite(
        network.config.energy_absolute_tolerance_j
            + network.config.energy_relative_tolerance * scale,
    )?;
    if residual.abs() > tolerance {
        return Err(ModalCouplingError::EnergyBalance {
            residual_j: residual,
            tolerance_j: tolerance,
        });
    }
    poll(gate)?;
    network.staged_frame = CoupledModalFrame {
        sample,
        observer_pressure_pa: pressure,
        modal_energy_j: modal_after,
        connection_energy_j: spring_after,
        external_work_j: work,
        component_dissipation_j: component_loss,
        connection_dissipation_j: connection_loss,
        energy_residual_j: residual,
        energy_tolerance_j: tolerance,
        solve_relative_residual: 0.0,
        connection_forces_n: network.reactions.clone(),
    };
    // Independent port observations cannot turn an accepted physical step into
    // a refusal. Preserve their original diagnostic instead of a zero account.
    let transfers = match &network.exact {
        Some(exact) => exact.connection_work(
            &x,
            &u,
            &rec_x,
            network.dt,
            &network.columns,
            &network.connections,
            connection_loss,
            gate,
        ),
        // None is admitted only when every connection has K=C=0. These
        // columns exert identically zero force throughout the actual step.
        None => Ok(vec![
            crate::port_adoption::PortTransfer::work(0.0, false);
            network.columns.len()
        ]),
    }
    .map_err(|error| fs_ports::CoupleError::ModuleObservationFailed {
        owner: "coupled-modal".into(),
        reason: error.to_string(),
    });
    network.staged_ports.accept_result(
        crate::port_adoption::StepAccount {
            tick: sample,
            duration_s: network.dt,
            storage_change_j: Some(after - before),
            loss_j: Some(component_loss + connection_loss),
            included_work_j: work,
            entropy_change_j_k: 0.0,
            tolerance_j: Some(tolerance),
            unavailable: None,
        },
        transfers,
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use fs_math::c64::C64;
    use fs_ports::PortedModule;

    fn network(rate: u32, omega: f64) -> CoupledModalSystem {
        let models: Vec<_> = [(omega, 0.01), (1.7 * omega, 0.2)]
            .into_iter()
            .enumerate()
            .map(|(i, (w, zeta))| {
                let mut model = ModalAcousticTimeModel::try_new(
                    rate,
                    vec![ModalAcousticMode {
                        angular_frequency_rad_s: w,
                        damping_ratio: zeta,
                        pressure_per_modal_velocity: C64::ZERO,
                    }],
                    ModalAcousticTimeBudget::audible_reference(),
                )
                .unwrap();
                model
                    .restore_states(&[ModalAcousticState {
                        displacement_m_sqrt_kg: if i == 0 { 1e-8 } else { -2e-8 },
                        velocity_m_sqrt_kg_per_s: 0.0,
                    }])
                    .unwrap();
                model
            })
            .collect();
        CoupledModalSystem::new(
            models,
            vec![ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: vec![1.0],
                },
                right: ModalAttachment {
                    component: 1,
                    shapes: vec![1.0],
                },
                stiffness_n_m: 0.2 * omega * omega,
                damping_n_s_m: 0.07 * omega,
                rest_extension_m: 1e-9,
            }],
            ModalCouplingConfig {
                max_modes: 2,
                max_connections: 1,
                max_setup_terms: 64,
                nyquist_guard_fraction: 0.9,
                maximum_total_energy_j: 1e6,
                maximum_abs_pressure_pa: 1e6,
                maximum_abs_connection_force_n: 1e6,
                solve_relative_tolerance: 1e-11,
                energy_absolute_tolerance_j: 1e-12,
                energy_relative_tolerance: 1e-10,
            },
            &CancelGate::new(),
        )
        .unwrap()
    }

    #[test]
    fn many_physical_dashpots_keep_exact_motion_work_and_caller_budgets() {
        // Seventy-two distinct no-slip/contact patches are physical ports,
        // even though all their displacements project onto one retained mode.
        // Their known aggregate gives an independent scalar ODE reference.
        let count = 72;
        let omega: f64 = 200.0;
        let q0: f64 = 0.001;
        let v0: f64 = 0.3;
        let damping_per_patch: f64 = 0.5;
        let mut model = ModalAcousticTimeModel::try_new(48_000,
            vec![ModalAcousticMode { angular_frequency_rad_s: omega,
                damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO }],
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q0,
            velocity_m_sqrt_kg_per_s: v0 }]).unwrap();
        let links: Vec<_> = (0..count).map(|_| ModalConnection {
            left: ModalAttachment { component: 0, shapes: vec![1.0] },
            right: ModalAttachment { component: 0, shapes: vec![0.0] },
            stiffness_n_m: 0.0, damping_n_s_m: damping_per_patch,
            rest_extension_m: 0.0,
        }).collect();
        let needed = (count + 1) * (count + 1);
        let config = ModalCouplingConfig {
            max_modes: 1, max_connections: count, max_setup_terms: needed,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
        };
        let gate = CancelGate::new_clock_free();
        // The cap generalization does not authorize extra caller work or ports.
        assert!(CoupledModalSystem::new(vec![model.clone()], links.clone(),
            ModalCouplingConfig { max_connections: count - 1, ..config }, &gate).is_err());
        assert!(CoupledModalSystem::new(vec![model.clone()], links.clone(),
            ModalCouplingConfig { max_setup_terms: needed - 1, ..config }, &gate).is_err());
        assert!(CoupledModalSystem::new(vec![model.clone()], links.clone(),
            ModalCouplingConfig { max_connections: MAX_MODAL_CONNECTIONS + 1, ..config }, &gate).is_err());
        let mut system = CoupledModalSystem::new(vec![model], links, config, &gate).unwrap();
        let gamma = 0.5 * count as f64 * damping_per_patch;
        let wd = (omega * omega - gamma * gamma).sqrt();
        let mut loss = 0.0;
        for step in 1..=8 {
            let frame = system.step(&[0.0]).unwrap().clone();
            let time = f64::from(step) / 48_000.0;
            let envelope = (-gamma * time).exp();
            let phase = wd * time;
            let q = envelope * (q0 * phase.cos() + (v0 + gamma * q0) / wd * phase.sin());
            let v = envelope * (v0 * phase.cos()
                - (gamma * v0 + omega * omega * q0) / wd * phase.sin());
            let state = system.components()[0].states()[0];
            assert!((state.displacement_m_sqrt_kg - q).abs() < 1e-12 * q.abs());
            assert!((state.velocity_m_sqrt_kg_per_s - v).abs() < 1e-12 * v.abs());
            assert_eq!(frame.connection_forces_n.len(), count);
            for &reaction in &frame.connection_forces_n {
                assert!((reaction + damping_per_patch * v).abs() < 1e-12 * v.abs());
            }
            let transfers = system.last_window_transfers().unwrap();
            assert_eq!(transfers.len(), count);
            let port_loss: f64 = transfers.iter().map(|t| -t.energy_into().value()).sum();
            assert!((port_loss - frame.connection_dissipation_j).abs()
                < 1e-11 * frame.connection_dissipation_j);
            assert!(frame.connection_dissipation_j > 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            loss += frame.connection_dissipation_j;
        }
        let initial = 0.5 * (omega * omega * q0 * q0 + v0 * v0);
        assert!((system.total_energy_j().unwrap() + loss - initial).abs() < 1e-11 * initial);
    }

    #[test]
    fn g0_disconnected_port_work_is_zero_on_the_accepted_factorized_step() {
        let seed = network(48_000, 1000.0);
        let mut network = CoupledModalSystem::new_with_columns(
            seed.models,
            vec![ModalColumnConnection {
                column: vec![1.0, -1.0],
                stiffness_n_m: 0.0,
                damping_n_s_m: 0.0,
                rest_extension_m: 1e-4,
            }],
            seed.config,
            &CancelGate::new(),
        )
        .unwrap();
        assert!(network.exact.is_none());
        assert!(network.last_window_transfers().is_none());
        network.step(&[1e-4, -2e-4]).unwrap();
        assert_eq!(network.samples_rendered(), 1);
        assert_eq!(network.last_window_duration_s(), Some(1.0 / 48_000.0));
        let work = network.last_window_transfers().unwrap();
        assert_eq!(work.len(), 1);
        assert_eq!(work[0].energy_into().value(), 0.0);
        assert!(network.last_window_transfer_error().is_none());
        let fs_ports::ModuleWindowAudit::Available(audit) = network.audit_last_window().unwrap()
        else {
            panic!("accepted disconnected window missing")
        };
        assert!(audit.violations().is_empty());
    }

    #[test]
    fn g3_connection_port_work_matches_exact_time_partition_and_body_balance() {
        let make = |rate| {
            let mut seed = network(rate, 1000.0);
            for (model, (q, v)) in seed.models.iter_mut().zip([(0.001, 0.3), (-0.002, -0.2)]) {
                model
                    .restore_states(&[ModalAcousticState {
                        displacement_m_sqrt_kg: q,
                        velocity_m_sqrt_kg_per_s: v,
                    }])
                    .unwrap();
            }
            CoupledModalSystem::new_with_columns(
                seed.models,
                vec![
                    ModalColumnConnection {
                        column: vec![1.0, -1.0],
                        stiffness_n_m: 2e5,
                        damping_n_s_m: 70.0,
                        rest_extension_m: 1e-4,
                    },
                    ModalColumnConnection {
                        column: vec![0.5, 1.0],
                        stiffness_n_m: 3e5,
                        damping_n_s_m: 25.0,
                        rest_extension_m: -2e-4,
                    },
                ],
                ModalCouplingConfig {
                    max_connections: 2,
                    max_setup_terms: 128,
                    ..seed.config
                },
                &CancelGate::new(),
            )
            .unwrap()
        };
        let mut whole = make(48_000);
        let mut half = make(96_000);
        let modal_before = component_energy(whole.components()).unwrap();
        let frame = whole.step(&[10.0, -20.0]).unwrap().clone();
        let work: Vec<_> = whole
            .last_window_transfers()
            .unwrap()
            .iter()
            .map(|t| t.energy_into().value())
            .collect();
        half.step(&[10.0, -20.0]).unwrap();
        let first = half.last_window_transfers().unwrap().to_vec();
        half.step(&[10.0, -20.0]).unwrap();
        for ((value, first), second) in work
            .iter()
            .zip(first)
            .zip(half.last_window_transfers().unwrap())
        {
            let split = first.energy_into().value() + second.energy_into().value();
            assert!((value - split).abs() <= 1e-9 * value.abs().max(split.abs()).max(1e-20));
        }
        let body_balance = frame.modal_energy_j - modal_before + frame.component_dissipation_j
            - frame.external_work_j;
        assert!((work.iter().sum::<f64>() - body_balance).abs() <= frame.energy_tolerance_j);
        let endpoint_product =
            frame.connection_forces_n[0] * dot(&whole.columns[0], &whole.free_delta).unwrap();
        assert!((work[0] - endpoint_product).abs() > 1e-8 * work[0].abs());
        let previous = whole.last_window_transfers().unwrap().to_vec();
        assert!(whole.step(&[f64::NAN, 0.0]).is_err());
        assert_eq!(whole.last_window_transfers().unwrap(), previous);
        let gate = CancelGate::new();
        gate.request();
        assert!(matches!(
            whole.step_under_gate(&[10.0, -20.0], &gate),
            Err(ModalCouplingError::Cancelled)
        ));
        assert_eq!(whole.last_window_transfers().unwrap(), previous);
    }

    #[test]
    fn complete_damped_pencil_preserves_time_partition_and_independent_losses() {
        // Unequal frequencies and losses make K and C noncommuting. Independent
        // component exponentials plus trapezoidal reactions fail this semigroup.
        for rate in [48_000, 96_000, 192_000] {
            for frequency in [1_000.0, 5_000.0, 10_000.0] {
                let omega = core::f64::consts::TAU * frequency;
                let mut whole = network(rate, omega);
                let mut half = network(2 * rate, omega);
                let frame = whole.step(&[1e-4, -2e-4]).unwrap().clone();
                half.step(&[1e-4, -2e-4]).unwrap();
                half.step(&[1e-4, -2e-4]).unwrap();
                let a = flatten(whole.components());
                let b = flatten(half.components());
                for (x, y) in a.iter().zip(b) {
                    assert!(
                        (x - y).abs() <= 1e-12 * x.abs().max(y.abs()).max(1e-12),
                        "{rate} {frequency}: {x} {y}"
                    );
                }
                assert!(
                    frame.component_dissipation_j > 0.0 && frame.connection_dissipation_j > 0.0
                );
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
    }

    #[test]
    fn full_network_capacity_and_cancel_refuse_without_dropping_state() {
        let model = ModalAcousticTimeModel::try_new(
            48_000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 1.0,
                damping_ratio: 0.01,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        let cap = (fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM - 2) / 6;
        let models = vec![model; cap + 1];
        let gate = CancelGate::new();
        // Capacity belongs to a genuinely connected exponential block.
        // Disconnected original coordinates must never consume that capacity.
        let links = [ConnectionLaw {
            stiffness_n_m: 1.0,
            damping_n_s_m: 0.0,
            rest_extension_m: 0.0,
        }];
        assert!(ExactNetwork::new(&models[..cap], &[vec![1.0; cap]], &links, cap, &gate).is_ok());
        assert!(matches!(
            ExactNetwork::new(&models, &[vec![1.0; cap + 1]], &links, cap + 1, &gate),
            Err(ModalCouplingError::Propagation(
                fs_phs::PhsError::LinearSplitBudget { .. }
            ))
        ));
        let independent = ExactNetwork::new(&models, &[], &[], cap + 1, &gate).unwrap();
        assert_eq!(independent.blocks.len(), cap + 1);
        let mut network = network(48_000, 1000.0);
        let before = flatten(network.components());
        gate.request();
        assert!(matches!(
            network.step_under_gate(&[0.0; 2], &gate),
            Err(ModalCouplingError::Cancelled)
        ));
        assert_eq!(flatten(network.components()), before);
        assert_eq!(network.samples_rendered(), 0);
    }

    #[test]
    fn noncontiguous_islands_match_the_complete_pencil_and_contact_response() {
        let base = network(48_000, 1_000.0);
        let models: Vec<_> = base.components().iter().cloned().cycle().take(4).collect();
        let modes: Vec<_> = models
            .iter()
            .flat_map(|m| m.modes().iter().copied())
            .collect();
        let columns = [vec![1.0, 0.0, -0.7, 0.0], vec![0.0, -0.3, 0.0, 1.0]];
        let links = [
            ConnectionLaw {
                stiffness_n_m: 2e5,
                damping_n_s_m: 70.0,
                rest_extension_m: 1e-9,
            },
            ConnectionLaw {
                stiffness_n_m: 3e5,
                damping_n_s_m: 90.0,
                rest_extension_m: -2e-9,
            },
        ];
        let gate = CancelGate::new();
        let islands = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
        assert_eq!(
            islands
                .blocks
                .iter()
                .map(|b| b.0.clone())
                .collect::<Vec<_>>(),
            vec![vec![0, 2], vec![1, 3]]
        );
        let full = ExactBlock::new(&modes, &columns, &links, &gate).unwrap();
        let dt = 1.0 / 48_000.0;
        let external = [1e-4, -2e-4, 3e-4, -4e-4];
        let (actual, x, u) = islands.advance(&models, &external, dt).unwrap();
        assert_eq!(u, full.held_input(&external).unwrap());
        let expected = step_render(&full.system, &full.split, &x, &u, dt).unwrap();
        // Energy changes, supplies and losses are differences of the stored
        // energy, so their roundoff scales with H rather than with themselves.
        let stored = full
            .system
            .hamiltonian(&x)
            .abs()
            .max(full.system.hamiltonian(&expected.x).abs());
        let close = |a: f64, b: f64| {
            assert!(
                (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(stored).max(1e-20),
                "{a:e} != {b:e} (stored energy {stored:e})"
            )
        };
        let roundoff = |a: f64, b: f64| {
            assert!(
                (a - b).abs() <= 256.0 * f64::EPSILON * a.abs().max(b.abs()).max(1e-20),
                "block/monolithic propagation beyond roundoff: {a:e} != {b:e}"
            )
        };
        for (&a, &b) in actual.x.iter().zip(&expected.x) {
            roundoff(a, b);
        }
        close(actual.delta_h, expected.delta_h);
        close(actual.supplied, expected.supplied);
        close(actual.dissipated, expected.dissipated);
        let losses = islands.losses(&x, &u, dt, actual.dissipated).unwrap();
        let full_losses = full.losses(&x, &u, dt, expected.dissipated).unwrap();
        assert!(losses.0 > 0.0 && losses.1 > 0.0);
        close(losses.0, full_losses.0);
        close(losses.1, full_losses.1);
        let work = islands
            .connection_work(
                &x,
                &u,
                &actual.x,
                dt,
                &columns,
                &links,
                losses.1,
                Some(&gate),
            )
            .unwrap();
        let full_work = full
            .connection_work(
                &x,
                &u,
                &expected.x,
                dt,
                &columns,
                &links,
                full_losses.1,
                Some(&gate),
            )
            .unwrap();
        for (a, b) in work.iter().zip(full_work) {
            close(a.energy_j, b.energy_j);
        }
        let contact = [1.0, -2.0, 0.3, 0.4];
        for (a, b) in islands
            .response(&contact, dt)
            .unwrap()
            .into_iter()
            .zip(full.response(&contact, dt).unwrap())
        {
            roundoff(a, b);
        }
        // Plant one physical cross-link, once as stiffness and once as loss.
        // Each must merge the two noncontiguous groups and change the motion;
        // fixture labels or an incorrect K-only graph would miss this.
        for (stiffness_n_m, damping_n_s_m) in [(4e5, 0.0), (0.0, 500.0)] {
            let cross_columns = [
                columns[0].clone(),
                columns[1].clone(),
                vec![0.0, 1.0, 1.0, 0.0],
            ];
            let cross_links = [
                ConnectionLaw {
                    stiffness_n_m: links[0].stiffness_n_m,
                    damping_n_s_m: links[0].damping_n_s_m,
                    rest_extension_m: links[0].rest_extension_m,
                },
                ConnectionLaw {
                    stiffness_n_m: links[1].stiffness_n_m,
                    damping_n_s_m: links[1].damping_n_s_m,
                    rest_extension_m: links[1].rest_extension_m,
                },
                ConnectionLaw {
                    stiffness_n_m,
                    damping_n_s_m,
                    rest_extension_m: 0.0,
                },
            ];
            let merged =
                ExactNetwork::new(&models, &cross_columns, &cross_links, 4, &gate).unwrap();
            assert_eq!(merged.blocks.len(), 1);
            assert_eq!(merged.blocks[0].0, vec![0, 1, 2, 3]);
            let (changed, _, held) = merged.advance(&models, &external, dt).unwrap();
            assert!(
                actual
                    .x
                    .iter()
                    .zip(&changed.x)
                    .any(|(a, b)| (a - b).abs() > 1e-6 * a.abs().max(b.abs()).max(1e-20)),
                "one cross-term must change the original-coordinate response"
            );
            let monolithic = ExactBlock::new(&modes, &cross_columns, &cross_links, &gate).unwrap();
            let expected =
                step_render(&monolithic.system, &monolithic.split, &x, &held, dt).unwrap();
            for (&a, &b) in changed.x.iter().zip(&expected.x) {
                roundoff(a, b);
            }
        }
    }
}
