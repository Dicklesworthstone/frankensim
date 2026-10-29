//! One continuous original modal/Maxwell/DG generator, cached at the fixed clock.
use super::*;
use crate::dg_port::{DgJunctionTrajectory, DgPortLoad, DgTerminalTrajectory, sum_terminal_work};
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
    junctions: Vec<JunctionCache>,
    owner_loss_integrals: Vec<usize>,
    column: Vec<f64>,
    pressure_column: Vec<f64>,
    storage: Vec<f64>,
    rest_force: Vec<f64>,
    link_integrals: Vec<Option<usize>>,
}
struct TerminalCache {
    /// Complete modal/DG/rational rows: constrained terminal flow can depend
    /// directly on the original inlet volume flow and hence modal velocity.
    flow_columns: Vec<Vec<f64>>,
    constraints: Vec<Vec<f64>>,
    constraint_norms: Vec<f64>,
    constraint_pressure_norms: Vec<f64>,
    extra_scales: Vec<f64>,
    work_start: usize,
    linear_start: usize,
}
impl TerminalCache {
    fn check_constraints(
        &self,
        state: &[f64],
        duct_dimension: usize,
    ) -> Result<(), ModalCouplingError> {
        if self.constraints.is_empty() {
            return Ok(());
        }
        if state.len() != duct_dimension + self.extra_scales.len() {
            return Err(invalid(
                "zero-direct terminal lost original constrained coordinates",
            ));
        }
        let mut square = 0.0;
        for (i, &value) in state.iter().enumerate() {
            let scaled = finite(
                value
                    * if i < duct_dimension {
                        1.0
                    } else {
                        self.extra_scales[i - duct_dimension]
                    },
            )?;
            square = finite(square + scaled * scaled)?;
        }
        let norm = det::sqrt(square);
        for (row, &dual_norm) in self.constraints.iter().zip(&self.constraint_norms) {
            let residual = dot(row, state)?;
            let scale = finite(dual_norm * norm)?;
            if residual.abs() > 1e-10 * scale.max(f64::MIN_POSITIVE) {
                return Err(invalid(
                    "zero-direct terminal state violates original pressure continuity; inconsistent states are never projected",
                ));
            }
        }
        Ok(())
    }
}
/// Pull a physical pressure row back through the existing exact solver's
/// numerical state scaling. This is an error scale, not terminal storage.
fn terminal_dual_norm(
    row: &[f64],
    duct_dimension: usize,
    extra_scales: &[f64],
) -> Result<f64, ModalCouplingError> {
    let mut square = 0.0;
    for (i, &value) in row.iter().enumerate() {
        let scaled = finite(
            value
                / if i < duct_dimension {
                    1.0
                } else {
                    extra_scales[i - duct_dimension]
                },
        )?;
        square = finite(square + scaled * scaled)?;
    }
    finite(det::sqrt(square))
}
struct JunctionCache {
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
        let nd = load.branches().iter().try_fold(
            load.duct().continuous_operator().dimension(),
            |dimension, branch| {
                dimension
                    .checked_add(branch.duct().continuous_operator().dimension())
                    .ok_or_else(|| invalid("joint branch dimension overflow"))
            },
        )?;
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
            + 6 * load.branches().len()
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
            || load.branches().iter().any(|branch| {
                branch.state().step_count() != 0
                    || branch.propagator().dt().to_bits() != network.dt.to_bits()
            })
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
        let dg = load.prepare_continuous_network(budget.max_dense_bytes, gate)?;
        let pressure_column = dg.inlet_pressure_column.clone();
        let mut j = vec![0.0; nd * nd];
        let mut r = vec![0.0; nd * nd];
        let mut q = vec![0.0; nd * nd];
        let mut g = vec![0.0; nd * count];
        let a = &dg.generator;
        for i in 0..nd {
            q[i * nd + i] = 1.0;
            for k in 0..nd {
                j[i * nd + k] = 0.5 * (a[i * nd + k] - a[k * nd + i]);
                r[i * nd + k] = -0.5 * (a[i * nd + k] + a[k * nd + i]);
            }
            for k in 0..count {
                g[i * count + k] = dg.inlet_column[i] * column[k];
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
            loss_form(embed(
                &dg.owners[0].loss,
                dg.owners[0].dimension,
                nm + dg.owners[0].offset,
            )),
            vec![0.0; n * n],
        ];
        // Integrate every original owner's physical loss separately. The
        // projected pHS resistance is used only to prepare the constrained
        // dynamics; it must not replace the actual local wall heat reports.
        let mut owner_loss_integrals = vec![2];
        for owner in &dg.owners[1..] {
            owner_loss_integrals.push(forms.len());
            forms.push(loss_form(embed(
                &owner.loss,
                owner.dimension,
                nm + owner.offset,
            )));
        }
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
        let attachment = load.prepare_continuous_terminal(&dg, budget.max_dense_bytes, gate)?;
        let mut terminal = if let Some(attachment) = &attachment {
            let work_start = forms.len();
            let linear_start = linear_forms.len();
            let mut flow_columns = Vec::with_capacity(attachment.ports.len());
            let mut constraint_pressure_norms = Vec::with_capacity(attachment.ports.len());
            for i in 0..attachment.ports.len() {
                let mut flow = vec![0.0; n];
                flow[nm..].copy_from_slice(&attachment.flow_columns[i]);
                for (k, &b) in column.iter().enumerate() {
                    flow[2 * k + 1] = attachment.inlet_flow_factors[i] * b;
                }
                for source in [
                    &attachment.pressure_columns[i],
                    &attachment.load_pressure_columns[i],
                ] {
                    let mut pressure = vec![0.0; n];
                    pressure[nm..].copy_from_slice(source);
                    let mut form = vec![0.0; n * n];
                    for a in 0..n {
                        for b in 0..n {
                            form[a * n + b] = 0.5 * (pressure[a] * flow[b] + flow[a] * pressure[b]);
                        }
                    }
                    forms.push(form);
                    linear_forms.push(pressure);
                }
                linear_forms.push(flow.clone());
                flow_columns.push(flow);
                constraint_pressure_norms.push(if attachment.constraints.is_empty() {
                    0.0
                } else {
                    let difference: Vec<_> = attachment.pressure_columns[i]
                        .iter()
                        .zip(&attachment.load_pressure_columns[i])
                        .map(|(a, b)| a - b)
                        .collect();
                    terminal_dual_norm(&difference, nd, &attachment.extra_scales)?
                });
            }
            let cache = TerminalCache {
                flow_columns,
                constraints: attachment.constraints.clone(),
                constraint_norms: attachment
                    .constraints
                    .iter()
                    .map(|row| terminal_dual_norm(row, nd, &attachment.extra_scales))
                    .collect::<Result<_, _>>()?,
                constraint_pressure_norms,
                extra_scales: attachment.extra_scales.clone(),
                work_start,
                linear_start,
            };
            let mut original_state = load.continuous_energy_coordinates()?;
            if let Some(original) = load.terminal() {
                original_state.extend(original.continuous_operator().real_coordinates());
            }
            cache.check_constraints(&original_state, nd)?;
            Some(cache)
        } else {
            None
        };
        let mut junctions = Vec::with_capacity(dg.junctions.len());
        for junction in &dg.junctions {
            poll(Some(gate))?;
            let mut flow = vec![0.0; n];
            flow[nm..nb].copy_from_slice(&junction.flow_column);
            for (i, &b) in column.iter().enumerate() {
                flow[2 * i + 1] = junction.inlet_flow_factor * b;
            }
            if let Some(terminal) = &terminal {
                for (factor, terminal_flow) in junction
                    .terminal_flow_factors
                    .iter()
                    .zip(&terminal.flow_columns)
                {
                    for (value, &coefficient) in flow.iter_mut().zip(terminal_flow) {
                        *value += factor * coefficient;
                    }
                }
            }
            let work_start = forms.len();
            let linear_start = linear_forms.len();
            for source in [
                &junction.bore_pressure_column,
                &junction.branch_pressure_column,
            ] {
                let mut pressure = vec![0.0; n];
                pressure[nm..nb].copy_from_slice(source);
                let mut form = vec![0.0; n * n];
                for a in 0..n {
                    for b in 0..n {
                        form[a * n + b] = 0.5 * (pressure[a] * flow[b] + flow[a] * pressure[b]);
                    }
                }
                forms.push(form);
                linear_forms.push(pressure);
            }
            linear_forms.push(flow);
            junctions.push(JunctionCache {
                work_start,
                linear_start,
            });
        }
        if let Some(terminal) = &mut terminal {
            terminal.linear_start += forms.len();
        }
        for junction in &mut junctions {
            junction.linear_start += forms.len();
        }
        let flow = if let Some(attachment) = attachment {
            let (j, r, _) = system.structure();
            let difference: Vec<_> = j.iter().zip(r).map(|(j, r)| j - r).collect();
            let mut generator = embed(&product(&difference, base_storage, nb), nb, 0);
            let width = nd + nt;
            for i in 0..width {
                generator[(nm + i) * n + nm..(nm + i) * n + n]
                    .copy_from_slice(&attachment.generator[i * width..(i + 1) * width]);
                if !attachment.constraints.is_empty() {
                    // Differentiated pressure closure changes the inlet vector
                    // as well as A. Retain the original mechanical reaction
                    // pressure; only the acoustic/rational receiving rows change.
                    for (k, &b) in column.iter().enumerate() {
                        generator[(nm + i) * n + 2 * k + 1] = attachment.inlet_column[i] * b;
                    }
                }
            }
            PreparedQuadraticFlow::new_extended_controlled(
                &system,
                &split,
                &generator,
                &attachment.extra_scales,
                network.dt,
                &forms,
                &linear_forms,
                budget,
                || gate.is_requested(),
            )
        } else {
            PreparedQuadraticFlow::new_controlled(
                &system,
                &split,
                network.dt,
                &forms,
                &linear_forms,
                budget,
                || gate.is_requested(),
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
            junctions,
            owner_loss_integrals,
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
            || load
                .branches()
                .iter()
                .any(|branch| branch.state().step_count() != network.samples_rendered())
        {
            return Err(invalid(
                "joint DG requires finite complete forcing and unchanged accepted clock",
            )
            .into());
        }
        let mut x = flatten(&network.models);
        x.extend_from_slice(&network.material_memory);
        x.extend(load.continuous_energy_coordinates()?);
        if let Some(terminal) = load.terminal() {
            x.extend(terminal.continuous_operator().real_coordinates());
        }
        if let Some(terminal) = &self.terminal {
            terminal.check_constraints(&x[self.modal_dimension..], self.duct_dimension)?;
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
        if let Some(terminal) = &self.terminal {
            terminal
                .check_constraints(&result.state[self.modal_dimension..], self.duct_dimension)?;
        }
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
        let (endpoint_input, observations) = if let Some(terminal) = &self.terminal {
            let endpoint_input: Vec<_> = terminal
                .flow_columns
                .iter()
                .map(|column| dot(column, &result.state))
                .collect::<Result<_, _>>()?;
            let observations: Vec<_> = (0..terminal.flow_columns.len())
                .map(|i| DgTerminalTrajectory {
                    pressure_impulse_pa_s: values[terminal.linear_start + 3 * i],
                    load_pressure_impulse_pa_s: values[terminal.linear_start + 3 * i + 1],
                    constraint_pressure_scale_pa_s: network.dt
                        * terminal.constraint_pressure_norms[i]
                        * root_bound,
                    outward_volume_m3: values[terminal.linear_start + 3 * i + 2],
                    port_work_j: values[terminal.work_start + 2 * i],
                    load_work_j: values[terminal.work_start + 2 * i + 1],
                })
                .collect();
            (endpoint_input, observations)
        } else {
            (Vec::new(), Vec::new())
        };
        let junctions: Vec<_> = self
            .junctions
            .iter()
            .map(|junction| DgJunctionTrajectory {
                bore_pressure_impulse_pa_s: values[junction.linear_start],
                branch_pressure_impulse_pa_s: values[junction.linear_start + 1],
                volume_m3: values[junction.linear_start + 2],
                bore_work_j: values[junction.work_start],
                branch_work_j: values[junction.work_start + 1],
            })
            .collect();
        let owner_losses: Vec<_> = self
            .owner_loss_integrals
            .iter()
            .map(|&index| values[index])
            .collect();
        let duct_loss = finite(owner_losses.iter().sum())?;
        let acoustic = load.preview_continuous_network_trajectory(
            &result.state[nm..nb],
            &result.state[nb..],
            &endpoint_input,
            network.dt,
            pressure_impulse,
            inlet_work_j,
            &observations,
            &junctions,
            &owner_losses,
            network.config.energy_absolute_tolerance_j,
            network.config.energy_relative_tolerance,
        )?;
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
        let loss = finite(modal.total_dissipation_j() + duct_loss)?;
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
