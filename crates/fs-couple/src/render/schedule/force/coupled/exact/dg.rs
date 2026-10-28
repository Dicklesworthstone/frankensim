//! One continuous original modal/Maxwell/DG generator, cached at the fixed clock.
use super::*;
use crate::dg_port::{DgPortLoad, DgTerminalTrajectory, sum_terminal_work};
use crate::render::schedule::force::coupled::dg::{ModalDgError, ModalDgFrame};
/// Explicit cold dimension (states + held inputs + one), bytes and work. No defaults.
pub use fs_phs::linear_exact::DenseExponentialBudget as ExactDgBudget;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;

pub(crate) struct PreparedJointDg {
    flow: PreparedQuadraticFlow,
    modal_dimension: usize,
    duct_dimension: usize,
    inlet_pressure_integral: usize,
    terminal: Option<TerminalCache>,
    column: Vec<f64>,
    pressure_column: Vec<f64>,
    storage: Vec<f64>,
    rest_force: Vec<f64>,
    link_integrals: Vec<Option<usize>>,
}
struct TerminalCache {
    flow_columns: Vec<Vec<f64>>,
    work_start: usize,
    linear_start: usize,
}
fn product(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                out[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    out
}
fn quadratic(a: &[f64], x: &[f64]) -> Result<f64, ModalCouplingError> {
    let n = x.len();
    let mut sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            sum = finite(sum + x[i] * a[i * n + j] * x[j])?;
        }
    }
    Ok(sum)
}
impl PreparedJointDg {
    pub(crate) fn new(
        network: &CoupledModalSystem,
        load: &DgPortLoad,
        column: Vec<f64>,
        budget: ExactDgBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalDgError> {
        poll(Some(gate))?;
        let dg = load.continuous_ordinary_duct()?;
        let nd = dg.dimension();
        let count = network.mode_count();
        let nm = 2 * count + network.material_memory.len();
        let nb = nm
            .checked_add(nd)
            .ok_or_else(|| invalid("joint DG dimension overflow"))?;
        let nt = load
            .terminal()
            .map_or(0, |t| t.continuous_operator().dimension());
        let n = nb
            .checked_add(nt)
            .ok_or_else(|| invalid("joint rational dimension overflow"))?;
        let d = n
            .checked_add(count)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| invalid("joint DG augmented dimension overflow"))?;
        let forms = 5
            + 5 * load.terminal_ports().len()
            + network
                .connections
                .iter()
                .filter(|c| c.damping_n_s_m > 0.0)
                .count();
        let bytes = d
            .checked_mul(d)
            .and_then(|v| v.checked_mul(8))
            .and_then(|v| v.checked_mul(2 * forms + 48))
            .ok_or_else(|| invalid("joint DG byte admission overflow"))?;
        let work = (d as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(256))
            .and_then(|v| v.checked_mul(forms as u64 + 1))
            .ok_or_else(|| invalid("joint DG work admission overflow"))?;
        if d > 256
            || d > budget.max_dimension
            || bytes > budget.max_dense_bytes
            || work > budget.max_multiply_terms
        {
            return Err(
                invalid("complete modal and DG operator exceeds explicit cold resources").into(),
            );
        }
        if network.samples_rendered() != 0
            || load.state().step_count() != 0
            || network.dt.to_bits() != load.propagator().dt().to_bits()
            || network.radiation.is_some()
            || column.len() != count
            || column.iter().any(|b| !b.is_finite())
            || network
                .connections
                .iter()
                .any(|c| c.constant_loss_reference_omega_rad_s.is_some())
        {
            return Err(invalid("exact DG requires cold matching clocks and ordinary modal/Maxwell source without direct air or constant-loss connections").into());
        }
        let modes: Vec<_> = network
            .models
            .iter()
            .flat_map(|m| m.modes().iter().copied())
            .collect();
        let base = ExactBlock::new_with_relaxation(
            &modes,
            &network.columns,
            &network.connections,
            network.relaxation_branches.clone(),
            gate,
        )?;
        let pressure_column = dg
            .ports()
            .iter()
            .find(|(id, _)| *id == load.inlet())
            .ok_or_else(|| invalid("original DG inlet missing"))?
            .1
            .clone();
        let mut j = vec![0.0; nd * nd];
        let mut r = vec![0.0; nd * nd];
        let mut q = vec![0.0; nd * nd];
        let mut g = vec![0.0; nd * count];
        let a = dg.generator_per_s();
        for i in 0..nd {
            q[i * nd + i] = 1.0;
            for k in 0..nd {
                j[i * nd + k] = 0.5 * (a[i * nd + k] - a[k * nd + i]);
                r[i * nd + k] = -0.5 * (a[i * nd + k] + a[k * nd + i]);
            }
            for k in 0..count {
                g[i * count + k] = pressure_column[i] * column[k];
            }
        }
        let duct = PortHamiltonian::new(
            nd,
            count,
            j,
            r.clone(),
            g,
            Box::new(QuadraticStorage::new(q, nd).map_err(ModalCouplingError::Propagation)?),
        )
        .map_err(ModalCouplingError::Propagation)?;
        let system = base
            .system
            .with_parallel_load(duct, &(0..count).collect::<Vec<_>>())
            .map_err(ModalCouplingError::Propagation)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: nb,
                max_omega_dt: core::f64::consts::PI,
            },
        )
        .map_err(ModalCouplingError::Propagation)?;
        let base_storage = split.quadratic_matrix();
        let embed = |source: &[f64], width: usize, offset: usize| {
            let mut out = vec![0.0; n * n];
            if !source.is_empty() {
                for i in 0..width {
                    out[(offset + i) * n + offset..(offset + i) * n + offset + width]
                        .copy_from_slice(&source[i * width..(i + 1) * width]);
                }
            }
            out
        };
        let storage = embed(base_storage, nb, 0);
        let loss_form =
            |resistance: Vec<f64>| product(&product(&storage, &resistance, n), &storage, n);
        let mut forms = vec![
            loss_form(embed(&base.component_r, nm, 0)),
            loss_form(embed(&base.material_r, nm, 0)),
            loss_form(embed(&r, nd, nm)),
            vec![0.0; n * n],
        ];
        for i in 0..count {
            for k in 0..nd {
                let value = 0.5 * column[i] * pressure_column[k];
                forms[3][(2 * i + 1) * n + nm + k] = value;
                forms[3][(nm + k) * n + 2 * i + 1] = value;
            }
        }
        let mut link_integrals = Vec::with_capacity(network.connections.len());
        for (b, law) in network.columns.iter().zip(&network.connections) {
            poll(Some(gate))?;
            if law.damping_n_s_m == 0.0 {
                link_integrals.push(None);
                continue;
            }
            let mut c = vec![0.0; n * n];
            for i in 0..count {
                for k in 0..count {
                    c[(2 * i + 1) * n + 2 * k + 1] = law.damping_n_s_m * b[i] * b[k];
                }
            }
            link_integrals.push(Some(forms.len()));
            forms.push(c);
        }
        let mut pressure = vec![0.0; n];
        pressure[nm..nm + nd].copy_from_slice(&pressure_column);
        let mut linear_forms = vec![pressure];
        let attachment = load.prepare_continuous_terminal(budget.max_dense_bytes, gate)?;
        let (flow, terminal) = if let Some(attachment) = attachment {
            let (j, r, _) = system.structure();
            let difference: Vec<_> = j.iter().zip(r).map(|(j, r)| j - r).collect();
            let mut generator = embed(&product(&difference, base_storage, nb), nb, 0);
            let width = nd + nt;
            for i in 0..width {
                generator[(nm + i) * n + nm..(nm + i) * n + n]
                    .copy_from_slice(&attachment.generator[i * width..(i + 1) * width]);
            }
            let work_start = forms.len();
            for i in 0..attachment.ports.len() {
                let flow = &attachment.flow_columns[i];
                for pressure in [
                    &attachment.pressure_columns[i],
                    &attachment.load_pressure_columns[i],
                ] {
                    let mut form = vec![0.0; n * n];
                    for a in 0..width {
                        for b in 0..width {
                            form[(nm + a) * n + nm + b] =
                                0.5 * (pressure[a] * flow[b] + flow[a] * pressure[b]);
                        }
                    }
                    forms.push(form);
                }
                for source in [
                    &attachment.pressure_columns[i],
                    &attachment.load_pressure_columns[i],
                    &attachment.flow_columns[i],
                ] {
                    let mut observation = vec![0.0; n];
                    observation[nm..].copy_from_slice(source);
                    linear_forms.push(observation);
                }
            }
            let terminal = TerminalCache {
                flow_columns: attachment.flow_columns,
                work_start,
                linear_start: forms.len() + 1,
            };
            let flow = PreparedQuadraticFlow::new_extended_controlled(
                &system,
                &split,
                &generator,
                &attachment.extra_scales,
                network.dt,
                &forms,
                &linear_forms,
                budget,
                || gate.is_requested(),
            );
            (flow, Some(terminal))
        } else {
            (
                PreparedQuadraticFlow::new_controlled(
                    &system,
                    &split,
                    network.dt,
                    &forms,
                    &linear_forms,
                    budget,
                    || gate.is_requested(),
                ),
                None,
            )
        };
        let flow = flow.map_err(|error| match error {
            fs_phs::PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
            fs_phs::PreparedStepError::Solver(error) => ModalCouplingError::Propagation(error),
        })?;
        poll(Some(gate))?;
        limit(
            "modal and DG total energy",
            finite(network.total_energy_j()? + load.stored_energy_j())?,
            network.config.maximum_total_energy_j,
        )?;
        Ok(Self {
            flow,
            modal_dimension: nm,
            duct_dimension: nd,
            inlet_pressure_integral: forms.len(),
            terminal,
            column,
            pressure_column,
            storage,
            rest_force: base.rest_force,
            link_integrals,
        })
    }

    pub(crate) fn step(
        &self,
        network: &mut CoupledModalSystem,
        load: &mut DgPortLoad,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<ModalDgFrame, ModalDgError> {
        poll(Some(gate))?;
        if external.len() != self.column.len()
            || external.iter().any(|f| !f.is_finite())
            || network.samples_rendered() != load.state().step_count()
        {
            return Err(invalid(
                "joint DG requires finite complete forcing and unchanged accepted clock",
            )
            .into());
        }
        let mut x = flatten(&network.models);
        x.extend_from_slice(&network.material_memory);
        x.extend_from_slice(
            load.continuous_ordinary_duct()?
                .energy_coordinates(load.state())
                .map_err(|e| {
                    ModalDgError::Acoustic(crate::acoustic_realize::AcousticRealizeError::Duct(e))
                })?,
        );
        if let Some(terminal) = load.terminal() {
            x.extend(terminal.continuous_operator().real_coordinates());
        }
        let input: Vec<_> = external
            .iter()
            .zip(&self.rest_force)
            .map(|(f, r)| finite(f + r))
            .collect::<Result<_, _>>()?;
        let before = finite(network.total_energy_j()? + load.stored_energy_j())?;
        let h0 = finite(0.5 * quadratic(&self.storage, &x)?)?;
        let result = self
            .flow
            .evaluate(&x, &input)
            .map_err(ModalCouplingError::Propagation)?;
        poll(Some(gate))?;
        let values = &result.integrals;
        let nm = self.modal_dimension;
        let inlet_work_j = values[3];
        let pressure_impulse = values[self.inlet_pressure_integral];
        let mean_pressure_pa = finite(pressure_impulse / network.dt)?;
        let nb = nm + self.duct_dimension;
        let endpoint_pressure = dot(&self.pressure_column, &result.state[nm..nb])?;
        let mut authored_external_work_j = 0.0;
        let mut volume = 0.0;
        for i in 0..external.len() {
            let dq = finite(result.state[2 * i] - x[2 * i])?;
            authored_external_work_j = finite(authored_external_work_j + external[i] * dq)?;
            volume = finite(volume + self.column[i] * dq)?;
        }
        let mut connection_loss = 0.0;
        let mut transfers = Vec::with_capacity(network.connections.len());
        for (index, (b, law)) in network.columns.iter().zip(&network.connections).enumerate() {
            let loss = self.link_integrals[index].map_or(0.0, |i| values[i]);
            connection_loss = finite(connection_loss + loss)?;
            let mut a = -law.rest_extension_m;
            let mut z = a;
            for i in 0..b.len() {
                a = finite(a + b[i] * x[2 * i])?;
                z = finite(z + b[i] * result.state[2 * i])?;
            }
            transfers.push(crate::port_adoption::PortTransfer::work(
                finite(-0.5 * law.stiffness_n_m * (z - a) * (z + a) - loss)?,
                false,
            ));
        }
        let root_bound = if self.terminal.is_some() {
            self.flow
                .interval_coordinate_norm_bound(&x, &input)
                .map_err(ModalCouplingError::Propagation)?
        } else {
            finite(det::sqrt(2.0 * h0) + network.dt * det::sqrt(input.iter().map(|f| f * f).sum()))?
        };
        let acoustic = if let Some(terminal) = &self.terminal {
            let endpoint_input: Vec<_> = terminal
                .flow_columns
                .iter()
                .map(|column| dot(column, &result.state[nm..]))
                .collect::<Result<_, _>>()?;
            let observations: Vec<_> = (0..terminal.flow_columns.len())
                .map(|i| DgTerminalTrajectory {
                    pressure_impulse_pa_s: values[terminal.linear_start + 3 * i],
                    load_pressure_impulse_pa_s: values[terminal.linear_start + 3 * i + 1],
                    outward_volume_m3: values[terminal.linear_start + 3 * i + 2],
                    port_work_j: values[terminal.work_start + 2 * i],
                    load_work_j: values[terminal.work_start + 2 * i + 1],
                })
                .collect();
            load.preview_continuous_terminal_trajectory(
                result.state[nm..nb].to_vec(),
                &result.state[nb..],
                &endpoint_input,
                network.dt,
                pressure_impulse,
                inlet_work_j,
                &observations,
                values[2],
                network.config.energy_absolute_tolerance_j,
                network.config.energy_relative_tolerance,
            )?
        } else {
            load.preview_continuous_trajectory(
                result.state[nm..nb].to_vec(),
                network.dt,
                pressure_impulse,
                inlet_work_j,
                values[2],
                network.config.energy_absolute_tolerance_j,
                network.config.energy_relative_tolerance,
            )?
        };
        let pressure_bound =
            finite(det::sqrt(self.pressure_column.iter().map(|b| b * b).sum()) * root_bound)?;
        limit(
            "DG inlet pressure over complete exact interval",
            pressure_bound,
            network.config.maximum_abs_pressure_pa,
        )?;
        let endpoint_forces: Vec<_> = external
            .iter()
            .zip(&self.column)
            .map(|(f, b)| finite(f - b * endpoint_pressure))
            .collect::<Result<_, _>>()?;
        stage_candidate(
            network,
            &endpoint_forces,
            PropagatedCandidate {
                endpoint: result.state[..nm].to_vec(),
                start: x[..nm].to_vec(),
                input,
                losses: [values[0], connection_loss, values[1], 0.0],
                // Without terminal storage, supply a conservative bound for the
                // complete interval; the original additional drive term is safe.
                initial_bound_energy: if self.terminal.is_some() {
                    finite(0.5 * root_bound.powi(2))?
                } else {
                    h0
                },
                external_work: Some(finite(authored_external_work_j - inlet_work_j)?),
                connection_transfers: Some(transfers),
            },
            Some(gate),
        )?;
        let modal = network.staged_frame.clone();
        let stored_energy_j =
            finite(modal.total_energy_j() + acoustic.output().ledger.stored_energy_j)?;
        let storage_change_j = finite(stored_energy_j - before)?;
        let loss = finite(modal.total_dissipation_j() + values[2])?;
        let terminal_work_j = sum_terminal_work(acoustic.terminal_frames());
        let residual =
            finite(storage_change_j + loss + terminal_work_j - authored_external_work_j)?;
        let tolerance = finite(
            network.config.energy_absolute_tolerance_j
                + network.config.energy_relative_tolerance
                    * (before.abs()
                        + stored_energy_j.abs()
                        + authored_external_work_j.abs()
                        + loss.abs()
                        + terminal_work_j.abs()),
        )?;
        limit(
            "modal and DG total energy",
            stored_energy_j,
            network.config.maximum_total_energy_j,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let port_work_residual_j = finite(acoustic.inlet_work_j() - inlet_work_j)?;
        let mean_flow_m3_s = finite(volume / network.dt)?;
        poll(Some(gate))?;
        let duct = acoustic.accept();
        network.publish_staged();
        Ok(ModalDgFrame {
            modal,
            duct,
            mean_flow_m3_s,
            mean_pressure_pa,
            stored_energy_j,
            storage_change_j,
            authored_external_work_j,
            inferred_total_dissipation_j: loss,
            terminal_work_j,
            port_work_residual_j,
            balance_residual_j: residual,
            balance_tolerance_j: tolerance,
            inlet_work_j,
            exact_joint_evolution: true,
        })
    }
}
