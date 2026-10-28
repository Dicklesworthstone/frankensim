//! Continuous single-duct composition keeps the original source and transaction.
use super::*;
use fs_duct::dg::continuous::{DgContinuous, state::DgTrajectoryPort};

impl DgPortLoad {
    pub(crate) fn continuous_ordinary_duct(
        &self,
    ) -> Result<DgContinuous<'_>, AcousticRealizeError> {
        if !self.branches.is_empty()
            || !self.terminal_channels.is_empty()
            || self.has_driven_mouth()
        {
            return Err(invalid(
                "continuous DG requires one original duct with no branch, gap or driven mouth",
            ));
        }
        let original = self.duct.continuous_operator();
        if original.ports().len() != 1 + self.terminal_ports().len()
            || self
                .terminal
                .as_ref()
                .is_some_and(|t| t.ports() != self.terminal_ports().len())
        {
            return Err(invalid(
                "continuous DG must retain every original terminal port",
            ));
        }
        Ok(original)
    }

    pub(crate) fn continuous_single_duct(&self) -> Result<DgContinuous<'_>, AcousticRealizeError> {
        if self.terminal.is_some()
            || !self.branches.is_empty()
            || !self.terminal_channels.is_empty()
            || self.has_driven_mouth()
            || self.duct.continuous_operator().ports().len() != 1
        {
            return Err(invalid(
                "continuous DG composition requires one original duct and inlet, with no rational terminal, branch, gap or driven mouth",
            ));
        }
        Ok(self.duct.continuous_operator())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn preview_continuous_trajectory(
        &mut self,
        endpoint: Vec<f64>,
        dt: f64,
        pressure_impulse: f64,
        inlet_work_j: f64,
        dissipated_j: f64,
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        let original = self.continuous_single_duct()?;
        let (state, output) = original
            .candidate_from_trajectory(
                &self.state,
                endpoint,
                dt,
                &[DgTrajectoryPort {
                    id: self.inlet,
                    pressure_impulse_pa_s: pressure_impulse,
                    work_into_duct_j: inlet_work_j,
                }],
                dissipated_j,
                absolute_tolerance_j,
                relative_tolerance,
            )
            .map_err(AcousticRealizeError::Duct)?;
        Ok(DgPortTrial {
            owner: self,
            state,
            output,
            terminal: None,
            terminal_frame: None,
            terminal_frames: Vec::new(),
            inlet_work_j,
            source_work_j: inlet_work_j,
            driven_mouth_frame: None,
            branch_states: Vec::new(),
            branch_frames: Vec::new(),
            terminal_channels: Vec::new(),
            terminal_channel_frames: Vec::new(),
        })
    }
}

/// Original DG/rational ODE after eliminating an invertible instantaneous D.
/// Terminal coordinates carry no declared storage; port work stays external.
pub(crate) struct DgContinuousTerminal {
    pub generator: Vec<f64>,
    pub extra_scales: Vec<f64>,
    pub ports: Vec<PortId>,
    pub flow_columns: Vec<Vec<f64>>,
    pub pressure_columns: Vec<Vec<f64>>,
    pub load_pressure_columns: Vec<Vec<f64>>,
}
impl DgPortLoad {
    pub(crate) fn prepare_continuous_terminal(
        &self,
        maximum_dense_bytes: usize,
        gate: &fs_exec::CancelGate,
    ) -> Result<Option<DgContinuousTerminal>, AcousticRealizeError> {
        let original = self.continuous_ordinary_duct()?;
        let Some(terminal) = &self.terminal else {
            return Ok(None);
        };
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let view = terminal.continuous_operator();
        let (nd, nt, p) = (original.dimension(), view.dimension(), terminal.ports());
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
        let solve = super::terminal::TerminalSolve::new(&model.d, p).map_err(|_| {
            invalid(
                "continuous rational terminal needs nonsingular well-conditioned instantaneous D",
            )
        })?;
        let ports = self.terminal_ports().to_vec();
        let mut pressure_columns = Vec::with_capacity(p);
        for id in &ports {
            let (_, column) = original
                .ports()
                .iter()
                .find(|(port, _)| port == id)
                .ok_or_else(|| invalid("continuous terminal lost its original DG port"))?;
            let mut pressure = vec![0.0; n];
            pressure[..nd].copy_from_slice(column);
            pressure_columns.push(pressure);
        }
        // Outward U = D^-1 (B_terminal^T x_DG - C z), in original order.
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
        let mut generator = vec![0.0; n * n];
        for i in 0..nd {
            generator[i * n..i * n + nd]
                .copy_from_slice(&original.generator_per_s()[i * nd..(i + 1) * nd]);
        }
        for i in 0..nt {
            generator[(nd + i) * n + nd..(nd + i) * n + n]
                .copy_from_slice(&model.a[i * nt..(i + 1) * nt]);
        }
        for i in 0..n {
            if gate.is_requested() {
                return Err(AcousticRealizeError::Cancelled);
            }
            for j in 0..n {
                for k in 0..p {
                    let input = if i < nd {
                        -pressure_columns[k][i]
                    } else {
                        model.b[(i - nd) * p + k]
                    };
                    generator[i * n + j] += input * flow_columns[k][j];
                }
            }
        }
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
    pub outward_volume_m3: f64,
    pub port_work_j: f64,
    pub load_work_j: f64,
}
impl DgPortLoad {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn preview_continuous_terminal_trajectory(
        &mut self,
        endpoint: Vec<f64>,
        terminal_endpoint: &[f64],
        terminal_endpoint_input: &[f64],
        dt: f64,
        pressure_impulse: f64,
        inlet_work_j: f64,
        terminals: &[DgTerminalTrajectory],
        dissipated_j: f64,
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        let original = self.continuous_ordinary_duct()?;
        if !dt.is_finite()
            || dt.to_bits() != self.propagator.dt().to_bits()
            || terminals.len() != self.terminal_ports().len()
            || !absolute_tolerance_j.is_finite()
            || absolute_tolerance_j < 0.0
            || !relative_tolerance.is_finite()
            || relative_tolerance < 0.0
        {
            return Err(invalid(
                "continuous terminal trajectory must retain clock, ports and finite tolerances",
            ));
        }
        let terminal = self
            .terminal
            .as_ref()
            .ok_or_else(|| invalid("continuous terminal owner missing"))?
            .continuous_operator()
            .candidate_from_coordinates(terminal_endpoint, terminal_endpoint_input)
            .map_err(|error| invalid(error.0))?;
        let mut ports = vec![DgTrajectoryPort {
            id: self.inlet,
            pressure_impulse_pa_s: pressure_impulse,
            work_into_duct_j: inlet_work_j,
        }];
        let mut terminal_frames = Vec::with_capacity(terminals.len());
        for (&id, data) in self.terminal_ports().iter().zip(terminals) {
            let residual = data.port_work_j - data.load_work_j;
            let work_scale = data.port_work_j.abs() + data.load_work_j.abs();
            let pressure_scale =
                data.pressure_impulse_pa_s.abs() + data.load_pressure_impulse_pa_s.abs();
            let flow = data.outward_volume_m3 / dt;
            let pressure = data.pressure_impulse_pa_s / dt;
            if ![residual, work_scale, pressure_scale, flow, pressure]
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
            ports.push(DgTrajectoryPort {
                id,
                pressure_impulse_pa_s: data.pressure_impulse_pa_s,
                work_into_duct_j: -data.port_work_j,
            });
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
        let (state, output) = original
            .candidate_from_trajectory(
                &self.state,
                endpoint,
                dt,
                &ports,
                dissipated_j,
                absolute_tolerance_j,
                relative_tolerance,
            )
            .map_err(AcousticRealizeError::Duct)?;
        let terminal_frame = match self.duct.outlet() {
            DgBoundary::Port(id) => terminal_frames.iter().find(|f| f.port == id).copied(),
            _ => None,
        };
        Ok(DgPortTrial {
            owner: self,
            state,
            output,
            terminal: Some(terminal),
            terminal_frame,
            terminal_frames,
            inlet_work_j,
            source_work_j: inlet_work_j,
            driven_mouth_frame: None,
            branch_states: Vec::new(),
            branch_frames: Vec::new(),
            terminal_channels: Vec::new(),
            terminal_channel_frames: Vec::new(),
        })
    }
}
