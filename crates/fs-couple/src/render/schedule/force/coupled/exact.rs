//! Full retained bilateral pencil; nonlinear contacts keep their existing
//! held discrete-gradient reactions around this SAME affine response.
use super::*;
use fs_math::det;
use fs_phs::linear_exact::{
    LinearSplit, LinearSplitBudget, integrated_resistance_loss, split_quadratic, step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, StepRecord};

pub(super) struct ExactNetwork {
    system: PortHamiltonian,
    split: LinearSplit,
    rest_force: Vec<f64>,
    component_r: Vec<f64>,
    connection_r: Vec<f64>,
}

impl ExactNetwork {
    pub(super) fn new(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        count: usize,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
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
        for model in models {
            poll(Some(gate))?;
            for mode in model.modes() {
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

    pub(super) fn advance(
        &self,
        models: &[ModalAcousticTimeModel],
        external: &[f64],
        dt: f64,
    ) -> Result<(StepRecord, Vec<f64>, Vec<f64>), ModalCouplingError> {
        let x = flatten(models);
        let u = self.held_input(external)?;
        let record = step_render(&self.system, &self.split, &x, &u, dt)
            .map_err(ModalCouplingError::Propagation)?;
        Ok((record, x, u))
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
        let h_x = exact.system.hamiltonian(&x);
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use fs_math::c64::C64;

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
        assert!(ExactNetwork::new(&models[..cap], &[], &[], cap, &gate).is_ok());
        assert!(matches!(
            ExactNetwork::new(&models, &[], &[], cap + 1, &gate),
            Err(ModalCouplingError::Propagation(
                fs_phs::PhsError::LinearSplitBudget { .. }
            ))
        ));
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
}
