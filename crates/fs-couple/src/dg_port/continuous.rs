//! Continuous original duct networks keep every source and transaction.
use super::*;
use fs_duct::dg::continuous::{DgContinuous, state::DgTrajectoryPort};

mod network;
mod zero_direct;
pub(crate) use network::DgContinuousNetwork;

/// Original DG/rational ODE after eliminating terminal flow multipliers.
/// Invertible D uses direct pressure closure; exactly zero D retains pressure
/// constraints and differentiates them once. Terminal coordinates carry no
/// declared storage; signed port work stays with the original external owner.
pub(crate) struct DgContinuousTerminal {
    pub generator: Vec<f64>,
    pub inlet_column: Vec<f64>,
    pub inlet_flow_factors: Vec<f64>,
    pub constraints: Vec<Vec<f64>>,
    pub extra_scales: Vec<f64>,
    pub ports: Vec<PortId>,
    pub flow_columns: Vec<Vec<f64>>,
    pub pressure_columns: Vec<Vec<f64>>,
    pub load_pressure_columns: Vec<Vec<f64>>,
}
impl DgPortLoad {
    pub(crate) fn prepare_continuous_terminal(
        &self,
        original: &DgContinuousNetwork,
        maximum_dense_bytes: usize,
        gate: &fs_exec::CancelGate,
    ) -> Result<Option<DgContinuousTerminal>, AcousticRealizeError> {
        let Some(terminal) = &self.terminal else {
            return Ok(None);
        };
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let view = terminal.continuous_operator();
        let (nd, nt, p) = (original.dimension, view.dimension(), terminal.ports());
        let n = nd
            .checked_add(nt)
            .ok_or_else(|| invalid("continuous terminal dimension overflow"))?;
        let bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(128))
            .ok_or_else(|| invalid("continuous terminal dense payload overflow"))?;
        if n > 256 || bytes > maximum_dense_bytes {
            return Err(invalid(
                "original DG and complete rational state exceed dense allowance",
            ));
        }
        let model = view
            .state_space(maximum_dense_bytes, &mut || gate.is_requested())
            .map_err(|error| {
                if gate.is_requested() {
                    AcousticRealizeError::Cancelled
                } else {
                    invalid(error.0)
                }
            })?;
        let zero_direct = model.d.iter().all(|d| *d == 0.0);
        let solve = if zero_direct {
            None
        } else {
            Some(super::terminal::TerminalSolve::new(&model.d, p).map_err(|_| {
                invalid(
                    "continuous rational terminal needs invertible D or exactly zero D; mixed-rank direct constraints are unsupported",
                )
            })?)
        };
        let ports = self.terminal_ports().to_vec();
        let mut pressure_columns = Vec::with_capacity(p);
        for column in &original.terminal_pressure_columns {
            let mut pressure = vec![0.0; n];
            pressure[..nd].copy_from_slice(column);
            pressure_columns.push(pressure);
        }
        let mut generator = vec![0.0; n * n];
        for i in 0..nd {
            generator[i * n..i * n + nd].copy_from_slice(&original.generator[i * nd..(i + 1) * nd]);
        }
        for i in 0..nt {
            generator[(nd + i) * n + nd..(nd + i) * n + n]
                .copy_from_slice(&model.a[i * nt..(i + 1) * nt]);
        }
        let mut inlet_column = vec![0.0; n];
        inlet_column[..nd].copy_from_slice(&original.inlet_column);
        let (flow_columns, inlet_flow_factors, constraints) = if let Some(solve) = solve {
            // Preserve direct closure and its original summation order:
            // outward U = D^-1 (B_terminal^T x_DG - C z).
            let mut flow_columns = vec![vec![0.0; n]; p];
            for col in 0..n {
                if gate.is_requested() {
                    return Err(AcousticRealizeError::Cancelled);
                }
                let rhs: Vec<_> = (0..p)
                    .map(|row| {
                        if col < nd {
                            pressure_columns[row][col]
                        } else {
                            -model.c[row * nt + col - nd]
                        }
                    })
                    .collect();
                for (row, value) in solve.apply(&rhs)?.into_iter().enumerate() {
                    flow_columns[row][col] = value;
                }
            }
            for i in 0..n {
                if gate.is_requested() {
                    return Err(AcousticRealizeError::Cancelled);
                }
                for j in 0..n {
                    for k in 0..p {
                        let input = if i < nd {
                            -original.terminal_input_columns[k][i]
                        } else {
                            model.b[(i - nd) * p + k]
                        };
                        generator[i * n + j] += input * flow_columns[k][j];
                    }
                }
            }
            (flow_columns, vec![0.0; p], Vec::new())
        } else {
            let mut inputs = vec![vec![0.0; n]; p];
            let mut constraints = pressure_columns.clone();
            for k in 0..p {
                for i in 0..nd {
                    inputs[k][i] = -original.terminal_input_columns[k][i];
                }
                for i in 0..nt {
                    inputs[k][nd + i] = model.b[i * p + k];
                    constraints[k][nd + i] = -model.c[k * nt + i];
                }
            }
            let retained =
                zero_direct::eliminate(&generator, &inputs, &constraints, &inlet_column, gate)?;
            generator = retained.generator;
            inlet_column = retained.inlet_column;
            (
                retained.flow_columns,
                retained.inlet_flow_factors,
                retained.constraints,
            )
        };
        let mut load_pressure_columns = vec![vec![0.0; n]; p];
        for i in 0..p {
            for j in 0..n {
                let mut value = if j < nd {
                    0.0
                } else {
                    model.c[i * nt + j - nd]
                };
                for k in 0..p {
                    value += model.d[i * p + k] * flow_columns[k][j];
                }
                load_pressure_columns[i][j] = value;
            }
        }
        // Shared scale for each original pole block balances its input/output
        // coefficients. It is numerical conditioning, never a stored-energy fit.
        let mut extra_scales = Vec::with_capacity(nt);
        let mut at = 0;
        for term in &terminal.source_model().poles {
            let width = term.state_dim() * p;
            let mut norm = 0.0_f64;
            for row in 0..p {
                for col in at..at + width {
                    norm = norm.max(model.c[row * nt + col].abs());
                }
            }
            let scale = if norm == 0.0 {
                1.0
            } else {
                fs_math::det::sqrt(norm)
            };
            extra_scales.extend(core::iter::repeat_n(scale, width));
            at += width;
        }
        if generator
            .iter()
            .chain(&inlet_column)
            .chain(&inlet_flow_factors)
            .chain(flow_columns.iter().flatten())
            .chain(load_pressure_columns.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "continuous terminal elimination left its finite source domain",
            ));
        }
        Ok(Some(DgContinuousTerminal {
            generator,
            inlet_column,
            inlet_flow_factors,
            constraints,
            extra_scales,
            ports,
            flow_columns,
            pressure_columns,
            load_pressure_columns,
        }))
    }
}

/// Complete independently integrated observations of one original terminal.
pub(crate) struct DgTerminalTrajectory {
    pub pressure_impulse_pa_s: f64,
    pub load_pressure_impulse_pa_s: f64,
    /// Original pressure-row norm times the exact interval coordinate bound;
    /// zero for direct closure. Keeps a zero-pressure constraint check scaled
    /// by physical operators when both measured impulses vanish.
    pub constraint_pressure_scale_pa_s: f64,
    pub outward_volume_m3: f64,
    pub port_work_j: f64,
    pub load_work_j: f64,
}

/// Independently integrated observations of both original junction traces.
pub(crate) struct DgJunctionTrajectory {
    pub bore_pressure_impulse_pa_s: f64,
    pub branch_pressure_impulse_pa_s: f64,
    pub volume_m3: f64,
    pub bore_work_j: f64,
    pub branch_work_j: f64,
}
impl DgPortLoad {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn preview_continuous_network_trajectory(
        &mut self,
        endpoint: &[f64],
        terminal_endpoint: &[f64],
        terminal_endpoint_input: &[f64],
        dt: f64,
        pressure_impulse: f64,
        inlet_work_j: f64,
        terminals: &[DgTerminalTrajectory],
        junctions: &[DgJunctionTrajectory],
        owner_losses: &[f64],
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        let dimension = self.admit_continuous_network()?;
        if endpoint.len() != dimension
            || endpoint.iter().any(|v| !v.is_finite())
            || !dt.is_finite()
            || dt.to_bits() != self.propagator.dt().to_bits()
            || terminals.len() != self.terminal_ports().len()
            || junctions.len() != self.branches.len()
            || owner_losses.len() != self.branches.len() + 1
            || !absolute_tolerance_j.is_finite()
            || absolute_tolerance_j < 0.0
            || !relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&relative_tolerance)
        {
            return Err(invalid(
                "continuous network trajectory must retain every original state, clock, port and finite tolerance",
            ));
        }
        let terminal = match &self.terminal {
            Some(original) => Some(
                original
                    .continuous_operator()
                    .candidate_from_coordinates(terminal_endpoint, terminal_endpoint_input)
                    .map_err(|error| invalid(error.0))?,
            ),
            None if terminal_endpoint.is_empty() && terminal_endpoint_input.is_empty() => None,
            None => {
                return Err(invalid(
                    "continuous terminal coordinates have no original owner",
                ));
            }
        };
        let mut ports = vec![DgTrajectoryPort {
            id: self.inlet,
            pressure_impulse_pa_s: pressure_impulse,
            work_into_duct_j: inlet_work_j,
        }];
        let mut terminal_frames = Vec::with_capacity(terminals.len());
        for (&id, data) in self.terminal_ports().iter().zip(terminals) {
            let residual = data.port_work_j - data.load_work_j;
            let work_scale = data.port_work_j.abs() + data.load_work_j.abs();
            let pressure_scale = data.pressure_impulse_pa_s.abs()
                + data.load_pressure_impulse_pa_s.abs()
                + data.constraint_pressure_scale_pa_s;
            let flow = data.outward_volume_m3 / dt;
            let pressure = data.pressure_impulse_pa_s / dt;
            if data.constraint_pressure_scale_pa_s < 0.0
                || ![residual, work_scale, pressure_scale, flow, pressure]
                    .iter()
                    .all(|v| v.is_finite())
                || residual.abs() > absolute_tolerance_j + relative_tolerance * work_scale
                || (data.pressure_impulse_pa_s - data.load_pressure_impulse_pa_s).abs()
                    > 1e-10 * pressure_scale.max(1e-300)
            {
                return Err(invalid(
                    "continuous original DG and rational terminal disagree on pressure or integrated work",
                ));
            }
            if let Some(index) = self.branches.iter().position(|branch| branch.port() == id) {
                let junction = &junctions[index];
                ports.push(DgTrajectoryPort {
                    id,
                    pressure_impulse_pa_s: junction.bore_pressure_impulse_pa_s,
                    work_into_duct_j: -junction.bore_work_j,
                });
            } else {
                ports.push(DgTrajectoryPort {
                    id,
                    pressure_impulse_pa_s: data.pressure_impulse_pa_s,
                    work_into_duct_j: -data.port_work_j,
                });
            }
            terminal_frames.push(DgTerminalFrame {
                port: id,
                outward_flow_m3_s: flow,
                mean_pressure_pa: pressure,
                port_work_j: data.port_work_j,
                load_work_j: data.load_work_j,
                channel_work_j: 0.0,
                port_work_residual_j: residual,
            });
        }
        let original = self.duct.continuous_operator();
        let main_dimension = original.dimension();
        let (state, mut output) = original
            .candidate_from_trajectory(
                &self.state,
                endpoint[..main_dimension].to_vec(),
                dt,
                &ports,
                owner_losses[0],
                absolute_tolerance_j,
                relative_tolerance,
            )
            .map_err(AcousticRealizeError::Duct)?;

        let mut branch_states = Vec::with_capacity(self.branches.len());
        let mut branch_frames = Vec::with_capacity(self.branches.len());
        let mut offset = main_dimension;
        for (index, (branch, junction)) in self.branches.iter().zip(junctions).enumerate() {
            let (DgBoundary::Port(inlet), DgBoundary::Port(outlet)) =
                (branch.duct().inlet(), branch.duct().outlet())
            else {
                return Err(invalid(
                    "continuous branch lost its two original flow ports",
                ));
            };
            let terminal_index = self
                .terminal_ports()
                .iter()
                .position(|id| *id == branch.port())
                .ok_or_else(|| invalid("continuous branch has no original distal terminal"))?;
            let external = &terminals[terminal_index];
            let view = branch.duct().continuous_operator();
            let end = offset + view.dimension();
            let (branch_state, branch_output) = view
                .candidate_from_trajectory(
                    branch.state(),
                    endpoint[offset..end].to_vec(),
                    dt,
                    &[
                        DgTrajectoryPort {
                            id: inlet,
                            pressure_impulse_pa_s: junction.branch_pressure_impulse_pa_s,
                            work_into_duct_j: junction.branch_work_j,
                        },
                        DgTrajectoryPort {
                            id: outlet,
                            pressure_impulse_pa_s: external.pressure_impulse_pa_s,
                            work_into_duct_j: -external.port_work_j,
                        },
                    ],
                    owner_losses[index + 1],
                    absolute_tolerance_j,
                    relative_tolerance,
                )
                .map_err(AcousticRealizeError::Duct)?;
            let bore_port = *output
                .port_outputs
                .iter()
                .find(|port| port.id == branch.port())
                .ok_or_else(|| invalid("continuous bore lost its original junction observation"))?;
            let branch_inlet = *branch_output
                .port_outputs
                .iter()
                .find(|port| port.id == inlet)
                .ok_or_else(|| invalid("continuous branch lost its original inlet observation"))?;
            let bore_trace = &original
                .ports()
                .iter()
                .find(|(id, _)| *id == branch.port())
                .ok_or_else(|| invalid("continuous bore lost its original junction trace"))?
                .1;
            let branch_trace = &view
                .ports()
                .iter()
                .find(|(id, _)| *id == inlet)
                .ok_or_else(|| invalid("continuous branch lost its original inlet trace"))?
                .1;
            // Component magnitudes give a cancellation-aware physical pressure
            // scale at the endpoint. Mean pressure alone can be zero during a
            // perfectly valid oscillation, so it is not a safe scale by itself.
            let endpoint_pressure_scale = bore_trace
                .iter()
                .zip(&endpoint[..main_dimension])
                .chain(branch_trace.iter().zip(&endpoint[offset..end]))
                .map(|(b, x)| (b * x).abs())
                .sum::<f64>();
            let impulse_scale = junction.bore_pressure_impulse_pa_s.abs()
                + junction.branch_pressure_impulse_pa_s.abs()
                + dt * endpoint_pressure_scale;
            let residual = junction.bore_work_j - junction.branch_work_j;
            let work_scale = junction.bore_work_j.abs() + junction.branch_work_j.abs();
            let flow = junction.volume_m3 / dt;
            if ![
                endpoint_pressure_scale,
                impulse_scale,
                residual,
                work_scale,
                flow,
            ]
            .iter()
            .all(|value| value.is_finite())
                || branch_output.ledger.step != output.ledger.step
                || (bore_port.end_conjugate - branch_inlet.end_conjugate).abs()
                    > 1e-10 * endpoint_pressure_scale.max(1e-300)
                || (junction.bore_pressure_impulse_pa_s - junction.branch_pressure_impulse_pa_s)
                    .abs()
                    > 1e-10 * impulse_scale.max(1e-300)
                || residual.abs() > absolute_tolerance_j + relative_tolerance * work_scale
            {
                return Err(invalid(
                    "continuous junction violates original pressure continuity or integrated work",
                ));
            }
            let distal = *branch_output
                .port_outputs
                .iter()
                .find(|port| port.id == outlet)
                .ok_or_else(|| {
                    invalid("continuous branch omitted its original distal observation")
                })?;
            let external_port = output
                .port_outputs
                .iter_mut()
                .find(|port| port.id == branch.port())
                .ok_or_else(|| invalid("continuous branch lost its external channel identity"))?;
            *external_port = fs_duct::dg::PortOutput {
                id: branch.port(),
                ..distal
            };
            let ledger = &branch_output.ledger;
            output.ledger.stored_energy_j += ledger.stored_energy_j;
            output.ledger.port_work_step_j += ledger.port_work_step_j;
            output.ledger.cumulative_port_work_j += ledger.cumulative_port_work_j;
            output.ledger.cumulative_parameter_work_j += ledger.cumulative_parameter_work_j;
            output.ledger.cumulative_dissipated_j += ledger.cumulative_dissipated_j;
            output.ledger.energy_balance_residual_j += ledger.energy_balance_residual_j;
            output.ledger.balance_checked &= ledger.balance_checked;
            branch_frames.push(DgBranchFrame {
                port: branch.port(),
                junction_flow_m3_s: flow,
                bore_port,
                branch_inlet,
                bore_work_j: junction.bore_work_j,
                branch_work_j: junction.branch_work_j,
                junction_work_residual_j: residual,
                duct: branch_output,
            });
            branch_states.push(branch_state);
            offset = end;
        }
        if !self.branches.is_empty() {
            output.ledger.relative_energy_drift = output.ledger.energy_balance_residual_j.abs()
                / (output.ledger.cumulative_port_work_j.abs()
                    + output.ledger.cumulative_parameter_work_j.abs())
                .max(1e-300);
            branch::validate_candidate(&state, &output, original.ports().len())?;
        }
        let terminal_frame = match self.duct.outlet() {
            DgBoundary::Port(id) => terminal_frames.iter().find(|f| f.port == id).copied(),
            _ => None,
        };
        Ok(DgPortTrial {
            owner: self,
            state,
            output,
            terminal,
            terminal_frame,
            terminal_frames,
            inlet_work_j,
            source_work_j: inlet_work_j,
            driven_mouth_frame: None,
            branch_states,
            branch_frames,
            terminal_channels: Vec::new(),
            terminal_channel_frames: Vec::new(),
        })
    }
}
