//! Frequency differentiation of the SAME complete source-bound DG-load pencil.
//! Every branch, rational state, channel cell and algebraic reaction remains.
use super::*;
use fs_duct::dg::{DgHarmonicDerivative, DgHarmonicPortDerivative};

/// Original rational terminal fields differentiated per rad/s. Generic terminal
/// heat/storage remains unavailable; signed pressure/flow work is observed.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicTerminalDerivative {
    pub ports: Vec<PortId>,
    pub state: Vec<C64>,
    pub outward_flow_m3_s: Vec<C64>,
    pub pressure_pa: Vec<C64>,
    pub mean_input_power_w: Vec<f64>,
    pub input_power_w: f64,
}

/// Original finite branch and zero-storage junction, differentiated per rad/s.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicBranchDerivative {
    pub port: PortId,
    pub junction_flow_m3_s: C64,
    pub bore_pressure_pa: C64,
    pub inlet_pressure_pa: C64,
    pub bore_power_w: f64,
    pub branch_power_w: f64,
    pub junction_power_residual_w: f64,
    pub duct: DgHarmonicDerivative,
}

/// Every original finite passage cell and independent storage/heat per rad/s.
#[derive(Clone, Debug, PartialEq)]
pub struct TerminalChannelHarmonicDerivative {
    pub velocities_m_s: Vec<C64>,
    pub pressure_drop_pa: C64,
    pub flow_m3_s: C64,
    pub mean_stored_energy_j: f64,
    pub viscous_dissipation_w: f64,
    pub input_power_w: f64,
    pub power_residual_w: f64,
    pub reactive_power_residual_w: f64,
    pub solve_relative_residual: f64,
}

/// Original finite passage/closure interface quantities, differentiated per rad/s.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicChannelDerivative {
    pub port: PortId,
    pub duct_pressure_pa: C64,
    pub exterior_pressure_pa: C64,
    pub outward_flow_m3_s: C64,
    pub pressure_drop_pa: C64,
    pub pressure_drop_residual_pa: C64,
    pub flow_residual_m3_s: C64,
    pub input_power_w: f64,
    pub interface_power_residual_w: f64,
    pub channel: Option<TerminalChannelHarmonicDerivative>,
}

/// Full original loaded response and its analytic frequency tangent.
///
/// Except `response` and `modal.response`, every physical field is a derivative
/// per rad/s. `modal` reuses the complete mechanical/material/air observation
/// owner; its embedded response equals `response.modal` and includes the actual
/// opposing duct reaction. Signed heat/storage slopes are not negative losses.
#[derive(Clone, Debug, PartialEq)]
pub struct DgLoadedModalHarmonicFrequencyResponse {
    pub response: DgLoadedModalHarmonicResponse,
    pub modal: ModalHarmonicFrequencyResponse,
    pub duct: DgHarmonicDerivative,
    pub branches: Vec<DgHarmonicBranchDerivative>,
    pub channels: Vec<DgHarmonicChannelDerivative>,
    pub terminal: Option<DgHarmonicTerminalDerivative>,
    pub terminal_power_w: f64,
    pub inlet_flow_m3_s: C64,
    pub inlet_pressure_pa: C64,
    pub authored_external_power_w: f64,
    pub mean_stored_energy_j: f64,
    pub interface_power_residual_w: f64,
    pub power_residual_w: f64,
    /// Backward residual of D*x' = F' - D'*x for the complete retained system.
    pub solve_relative_residual: f64,
}

impl DgLoadedModalHarmonicFrequencyResponse {
    /// Derivative of known physical losses, excluding unknown rational-terminal
    /// internal heat. Signed terminal input work remains a separate exchange.
    pub fn total_dissipation_derivative_w(&self) -> f64 {
        let duct_loss = |duct: &DgHarmonicDerivative| {
            duct.wall_dissipation_w + duct.molecular_dissipation_w + duct.matched_absorption_w
        };
        self.modal.total_dissipation_derivative()
            + duct_loss(&self.duct)
            + self
                .branches
                .iter()
                .map(|branch| duct_loss(&branch.duct))
                .sum::<f64>()
            + self
                .channels
                .iter()
                .map(|channel| {
                    channel
                        .channel
                        .as_ref()
                        .map_or(0.0, |field| field.viscous_dissipation_w)
                })
                .sum::<f64>()
    }
}

fn power_derivative(pressure: C64, flow: C64, dp: C64, dq: C64) -> f64 {
    0.5 * (dp * flow.conj() + pressure * dq.conj()).re
}

fn port_derivative(
    duct: &DgHarmonicDerivative,
    id: PortId,
) -> Result<&DgHarmonicPortDerivative, ModalCouplingError> {
    duct.ports
        .iter()
        .find(|port| port.id == id)
        .ok_or_else(|| invalid("original DG tangent observation omitted a declared port"))
}

impl CoupledModalSystem {
    /// Analytic frequency tangent of the complete bare-duct loaded response.
    /// Shares one joint LU and retains the original material/air and DG memory.
    /// No isolated duct impedance is inverted, including at a bare duct pole.
    /// `force_derivative` is dF/domega; None holds authored force phasors fixed.
    /// The original inlet projection and all physical coefficients stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn harmonic_response_with_duct_frequency_derivative(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicFrequencyResponse, DgHarmonicCouplingError> {
        self.harmonic_frequency_bound(
            omega_rad_s,
            force,
            force_derivative,
            duct,
            inlet,
            flow_column,
            None,
            &[],
            &[],
            maximum_work,
            gate,
        )
    }

    /// Analytic tangent of the COMPLETE source-bound loaded system. Original
    /// DG gas/wall/molecular coordinates, finite branches, MIMO terminal states,
    /// finite passage cells, junctions and rigid closure reactions all remain.
    /// Strictly proper/singular-direct terminals do not require any inverse.
    /// The same a+s*b assembly supplies acoustic derivative actions, and a
    /// second triangular solve reuses the original joint LU. No dense derivative
    /// matrix, fitted pole response or finite difference is used.
    ///
    /// Physical pressure/flow, source heat and storage slopes are independently
    /// mapped in their original owners. Unknown rational-terminal heat/storage
    /// is not invented. Source bands and the primary state's local Mach limit
    /// remain enforced; derivative amplitudes are slopes, not new fluid speeds.
    /// Driven mouths still require retained source equations and explicitly refuse.
    ///
    /// With d complete unknowns, C connections, G anchored-loss links, M Maxwell
    /// arms, f finite channel cells/drops and c channel inputs, work is bounded
    /// by d²*(d+G+2*(C+M+12))+512*f+128*d*c. Full payload is bounded before
    /// assembly at MAX_HARMONIC_BATCH_RESPONSE_BYTES; no state is truncated.
    /// No accepted physical state, sample clock or ledger changes.
    pub fn harmonic_response_with_load_frequency_derivative(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        load: &DgPortLoad,
        flow_column: &[f64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicFrequencyResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        if load.has_driven_mouth() {
            return Err(invalid(
                "continuous modal/DG tangent requires driven-mouth source equations to be retained",
            )
            .into());
        }
        self.harmonic_frequency_bound(
            omega_rad_s,
            force,
            force_derivative,
            load.duct(),
            load.inlet(),
            flow_column,
            load.terminal().map(|owner| (owner, load.terminal_ports())),
            load.branches(),
            load.terminal_channels(),
            maximum_work,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn harmonic_frequency_bound(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        terminal: Option<(&PreparedMimoLoad, &[PortId])>,
        branches: &[DgBranchLoad],
        channels: &[DgTerminalChannelLoad],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicFrequencyResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        let n = self.mode_count();
        if force_derivative.is_some_and(|values| {
            values.len() != n || values.iter().any(|&value| !complex_finite(value))
        }) {
            return Err(
                invalid("loaded force derivative must retain every finite original mode").into(),
            );
        }
        let view = duct.continuous_operator();
        let branch_views: Vec<_> = branches
            .iter()
            .map(|branch| branch.duct().continuous_operator())
            .collect();
        let binding = operator::DgBinding::new(
            n,
            &view,
            inlet,
            flow_column,
            terminal,
            branches,
            &branch_views,
            channels,
            gate,
        )?;
        binding.admit_frequency(omega_rad_s, gate)?;
        let extra = binding.additional_dimension()?;
        let air = self
            .radiation
            .as_ref()
            .map_or(0, |owner| owner.model.poles.len());
        let dimension = n
            .checked_add(air)
            .and_then(|value| value.checked_add(extra))
            .ok_or_else(|| invalid("complete DG tangent dimension overflow"))?;
        let channel_work = binding
            .channel_dimension
            .checked_mul(512)
            .and_then(|value| {
                dimension
                    .checked_mul(channels.len())
                    .and_then(|count| count.checked_mul(128))
                    .and_then(|channel| value.checked_add(channel))
            })
            .ok_or_else(|| invalid("complete channel tangent work overflow"))?;
        let work = harmonic_work(
            dimension,
            self.connections.len(),
            self.connections
                .iter()
                .filter(|law| law.constant_loss_reference_omega_rad_s.is_some())
                .count(),
            self.relaxation_branches.len(),
            12,
            2,
        )?
        .checked_add(channel_work)
        .ok_or_else(|| invalid("complete DG tangent work overflow"))?;
        // Three modal payloads include the paired native modal response. The
        // linear full-state allowance covers DG physical fields/ports, terminal
        // vectors, channel cells, every response header and tangent scratch.
        let bytes = harmonic_response_bytes(
            n,
            self.connections.len(),
            self.relaxation_branches.len(),
            air,
            3,
        )?
        .checked_add(
            dimension
                .checked_mul(1536)
                .ok_or_else(|| invalid("complete DG tangent payload overflow"))?,
        )
        .ok_or_else(|| invalid("complete DG tangent payload overflow"))?;
        if dimension > MAX_HARMONIC_UNKNOWNS
            || work > maximum_work
            || bytes > MAX_HARMONIC_BATCH_RESPONSE_BYTES
        {
            return Err(invalid(
                "complete DG tangent exceeds declared work or response payload capacity",
            )
            .into());
        }
        let mut operator =
            self.harmonic_operator(omega_rad_s, force, extra, maximum_work - channel_work, gate)?;
        let offset = operator.original_dimension;
        let appended = binding.append(
            offset,
            flow_column,
            16 * dimension * dimension,
            &mut operator::AcEntries {
                values: &mut operator.values,
                dimension,
                s: C64::new(0.0, -omega_rad_s),
            },
            gate,
        )?;
        let factor = HarmonicFactor::new(&operator.values, dimension, gate)?;
        let (state, relative) = factor.solve(force, self.config.solve_relative_tolerance, gate)?;
        let response = self.observe_harmonic_dg_state(
            omega_rad_s,
            force,
            inlet,
            flow_column,
            branches,
            &binding,
            &appended,
            &state,
            relative,
            maximum_work,
            gate,
        )?;
        let mut rhs =
            self.harmonic_frequency_rhs(omega_rad_s, &state[..offset], force_derivative, gate)?;
        rhs.resize(dimension, C64::ZERO);
        binding.apply_frequency_derivative(
            offset,
            flow_column,
            &appended,
            &state,
            &mut rhs,
            gate,
        )?;
        let (derivative, residual) =
            factor.solve(&rhs, self.config.solve_relative_tolerance, gate)?;
        self.observe_dg_frequency_derivative(
            response,
            force,
            force_derivative,
            inlet,
            flow_column,
            branches,
            &binding,
            &appended,
            &state,
            &derivative,
            residual,
            maximum_work,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_dg_frequency_derivative(
        &self,
        response: DgLoadedModalHarmonicResponse,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        inlet: PortId,
        flow_column: &[f64],
        branches: &[DgBranchLoad],
        binding: &operator::DgBinding<'_>,
        appended: &operator::AppendedDg,
        state: &[C64],
        derivative: &[C64],
        relative: f64,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicFrequencyResponse, DgHarmonicCouplingError> {
        let n = self.mode_count();
        let dimension = state.len();
        let offset = dimension - binding.additional_dimension()?;
        let omega = response.modal.omega_rad_s;
        let s = C64::new(0.0, -omega);
        let flow_derivative = s * projection(flow_column, &derivative[..n])
            + C64::new(0.0, -1.0) * projection(flow_column, &state[..n]);
        let mut inflows = vec![(inlet, flow_derivative)];
        for (j, &id) in binding.terminal_ports.iter().enumerate() {
            let outward = match branches.iter().position(|branch| branch.port() == id) {
                Some(index) => derivative[appended.junction_offset + index],
                None => derivative[appended.flow_offset + j],
            };
            inflows.push((id, -outward));
        }
        for (index, branch) in branches.iter().enumerate() {
            if binding.branch_ports[index].1.is_none() {
                inflows.push((branch.port(), -derivative[appended.junction_offset + index]));
            }
        }
        let duct = binding
            .view
            .observe_harmonic_frequency_derivative(
                &response.duct,
                &derivative[offset..offset + binding.view.dimension()],
                &inflows,
                maximum_work,
                self.config.solve_relative_tolerance,
                &mut || gate.is_requested(),
            )
            .map_err(|error| duct_error(error, gate))?;
        let pressure_derivative = port_derivative(&duct, inlet)?.pressure_pa;
        let applied: Vec<_> = force
            .iter()
            .zip(flow_column)
            .map(|(&f, &b)| f - response.inlet_pressure_pa.scale(b))
            .collect();
        let applied_derivative: Vec<_> = flow_column
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                force_derivative.map_or(C64::ZERO, |values| values[i])
                    - pressure_derivative.scale(b)
            })
            .collect();
        let modal = self.observe_harmonic_frequency_derivative(
            response.modal.clone(),
            &applied,
            Some(&applied_derivative),
            &state[..offset],
            &derivative[..offset],
            relative,
            gate,
        )?;
        let terminal = observe_terminal(
            binding,
            appended,
            &response,
            state,
            derivative,
            self.config.solve_relative_tolerance,
            gate,
        )?;
        let terminal_power_w = terminal
            .as_ref()
            .map_or(0.0, |terminal| terminal.input_power_w);
        let tolerance =
            self.config.solve_relative_tolerance + 64.0 * dimension as f64 * f64::EPSILON;
        let mut branch_derivatives = Vec::with_capacity(branches.len());
        for (index, (branch, view)) in branches.iter().zip(binding.branch_views).enumerate() {
            poll(Some(gate))?;
            let start = offset + binding.branch_offsets[index];
            let flow = response.branches[index].junction_flow_m3_s;
            let d_flow = derivative[appended.junction_offset + index];
            let mut inputs = vec![(binding.branch_ports[index].0, d_flow)];
            if let Some(outlet) = binding.branch_ports[index].1 {
                let port = binding
                    .terminal_ports
                    .iter()
                    .position(|id| *id == branch.port())
                    .ok_or_else(|| invalid("branch tangent lost its original outlet"))?;
                inputs.push((outlet, -derivative[appended.flow_offset + port]));
            }
            let branch_duct = view
                .observe_harmonic_frequency_derivative(
                    &response.branches[index].duct,
                    &derivative[start..start + view.dimension()],
                    &inputs,
                    maximum_work,
                    self.config.solve_relative_tolerance,
                    &mut || gate.is_requested(),
                )
                .map_err(|error| duct_error(error, gate))?;
            let bore_pressure = port_derivative(&duct, branch.port())?.pressure_pa;
            let inlet_pressure =
                port_derivative(&branch_duct, binding.branch_ports[index].0)?.pressure_pa;
            let original = &response.branches[index];
            let bore_power =
                power_derivative(original.bore_pressure_pa, flow, bore_pressure, d_flow);
            let branch_power =
                power_derivative(original.inlet_pressure_pa, flow, inlet_pressure, d_flow);
            let pressure_scale = |values: &[C64]| {
                binding.junction_columns[index]
                    .0
                    .iter()
                    .enumerate()
                    .map(|(i, b)| b.abs() * values[offset + i].abs())
                    .sum::<f64>()
                    + binding.junction_columns[index]
                        .1
                        .iter()
                        .enumerate()
                        .map(|(i, b)| b.abs() * values[start + i].abs())
                        .sum::<f64>()
            };
            let scale = pressure_scale(derivative);
            let power_scale = 0.5 * (scale * flow.abs() + pressure_scale(state) * d_flow.abs())
                + bore_power.abs()
                + branch_power.abs();
            if ![scale, power_scale, bore_power, branch_power]
                .iter()
                .all(|value| value.is_finite())
                || (bore_pressure - inlet_pressure).abs() > tolerance * scale.max(1e-300)
                || (bore_power - branch_power).abs() > tolerance * power_scale.max(1e-300)
            {
                return Err(
                    invalid("original branch tangent violates junction pressure or power").into(),
                );
            }
            branch_derivatives.push(DgHarmonicBranchDerivative {
                port: branch.port(),
                junction_flow_m3_s: d_flow,
                bore_pressure_pa: bore_pressure,
                inlet_pressure_pa: inlet_pressure,
                bore_power_w: bore_power,
                branch_power_w: branch_power,
                junction_power_residual_w: bore_power - branch_power,
                duct: branch_duct,
            });
        }
        let mut channel_derivatives = Vec::with_capacity(binding.channels.len());
        for (index, owner) in binding.channels.iter().enumerate() {
            poll(Some(gate))?;
            let port = owner.input().port;
            let port_index = binding
                .terminal_ports
                .iter()
                .position(|id| *id == port)
                .ok_or_else(|| invalid("channel tangent lost its original terminal row"))?;
            let exterior = terminal
                .as_ref()
                .ok_or_else(|| invalid("channel tangent requires the complete exterior"))?;
            let original = &response.channels[index];
            let (observed, observed_port) =
                match branches.iter().position(|branch| branch.port() == port) {
                    Some(branch) => (
                        &branch_derivatives[branch].duct,
                        binding.branch_ports[branch].1.ok_or_else(|| {
                            invalid("channel tangent cannot follow a sealed branch")
                        })?,
                    ),
                    None => (&duct, port),
                };
            let dp_duct = port_derivative(observed, observed_port)?.pressure_pa;
            let dp_exterior = exterior.pressure_pa[port_index];
            let d_drop = dp_duct - dp_exterior;
            let d_flow = exterior.outward_flow_m3_s[port_index];
            let drop = original.duct_pressure_pa - original.exterior_pressure_pa;
            let input_power = power_derivative(drop, original.outward_flow_m3_s, d_drop, d_flow);
            let (relative_offset, column) = binding.terminal_columns[port_index];
            let source = appended
                .terminal_source
                .as_ref()
                .ok_or_else(|| invalid("channel tangent lost source coefficients"))?;
            let pressure_scale = |values: &[C64]| {
                column
                    .iter()
                    .enumerate()
                    .map(|(i, b)| b.abs() * values[offset + relative_offset + i].abs())
                    .sum::<f64>()
                    + (0..binding.terminal_dimension)
                        .map(|i| {
                            source.c[port_index * binding.terminal_dimension + i].abs()
                                * values[appended.terminal_offset + i].abs()
                        })
                        .sum::<f64>()
                    + (0..binding.terminal_ports.len())
                        .map(|i| {
                            source.d[port_index * binding.terminal_ports.len() + i].abs()
                                * values[appended.flow_offset + i].abs()
                        })
                        .sum::<f64>()
            };
            let (constitutive, flow_residual, field, flow_scale, pressure_residual) =
                if let Some(channel) = owner.channel() {
                    let view = channel.continuous_operator()?;
                    let start = appended.channel_offsets[index];
                    let velocities = &derivative[start..start + view.dimension()];
                    let constitutive = derivative[start + view.dimension()];
                    let original_field = original
                        .channel
                        .as_ref()
                        .ok_or_else(|| invalid("channel tangent lost its original field"))?;
                    let field = observe_channel(
                        view,
                        original_field,
                        velocities,
                        d_drop,
                        omega,
                        tolerance,
                        gate,
                    )?;
                    let flow_scale = d_flow.abs()
                        + view
                            .cell_areas_m2()
                            .iter()
                            .zip(velocities)
                            .map(|(area, velocity)| area * velocity.abs())
                            .sum::<f64>();
                    (
                        constitutive,
                        field.flow_m3_s - d_flow,
                        Some(field),
                        flow_scale,
                        d_drop - constitutive,
                    )
                } else {
                    (
                        d_drop,
                        -d_flow,
                        None,
                        flow_derivative.abs()
                            + exterior
                                .outward_flow_m3_s
                                .iter()
                                .map(|value| value.abs())
                                .sum::<f64>(),
                        C64::ZERO,
                    )
                };
            let scale = pressure_scale(derivative) + constitutive.abs();
            let field_power = field.as_ref().map_or(0.0, |field| field.input_power_w);
            let work_scale = input_power.abs()
                + field_power.abs()
                + 0.5
                    * (scale * original.outward_flow_m3_s.abs()
                        + pressure_scale(state) * flow_scale);
            let work_residual = input_power - field_power;
            if ![
                dp_duct,
                dp_exterior,
                d_flow,
                constitutive,
                pressure_residual,
                flow_residual,
            ]
            .iter()
            .all(|&value| complex_finite(value))
                || ![scale, flow_scale, work_scale, work_residual]
                    .iter()
                    .all(|value| value.is_finite())
                || pressure_residual.abs() > tolerance * scale.max(1e-300)
                || flow_residual.abs() > tolerance * flow_scale.max(1e-300)
                || work_residual.abs() > tolerance * work_scale.max(1e-300)
            {
                return Err(invalid(
                    "original channel tangent disagrees on pressure, flow or power",
                )
                .into());
            }
            channel_derivatives.push(DgHarmonicChannelDerivative {
                port,
                duct_pressure_pa: dp_duct,
                exterior_pressure_pa: dp_exterior,
                outward_flow_m3_s: d_flow,
                pressure_drop_pa: constitutive,
                pressure_drop_residual_pa: pressure_residual,
                flow_residual_m3_s: flow_residual,
                input_power_w: input_power,
                interface_power_residual_w: work_residual,
                channel: field,
            });
        }
        let authored_external_power_w = force
            .iter()
            .enumerate()
            .map(|(i, &f)| {
                let df = force_derivative.map_or(C64::ZERO, |values| values[i]);
                power_derivative(
                    f,
                    response.modal.velocity_m_sqrt_kg_per_s[i],
                    df,
                    modal.velocity_derivative[i],
                )
            })
            .sum::<f64>();
        let mean_stored_energy_j = modal.mean_stored_energy_derivative
            + duct.mean_stored_energy_j
            + branch_derivatives
                .iter()
                .map(|branch| branch.duct.mean_stored_energy_j)
                .sum::<f64>()
            + channel_derivatives
                .iter()
                .map(|channel| {
                    channel
                        .channel
                        .as_ref()
                        .map_or(0.0, |field| field.mean_stored_energy_j)
                })
                .sum::<f64>();
        let interface_power_residual_w = authored_external_power_w
            - modal.external_power_derivative
            - duct.input_power_w
            - branch_derivatives
                .iter()
                .map(|branch| branch.duct.input_power_w)
                .sum::<f64>()
            - terminal_power_w
            - channel_derivatives
                .iter()
                .map(|channel| channel.input_power_w)
                .sum::<f64>();
        let mut result = DgLoadedModalHarmonicFrequencyResponse {
            response,
            modal,
            duct,
            branches: branch_derivatives,
            channels: channel_derivatives,
            terminal,
            terminal_power_w,
            inlet_flow_m3_s: flow_derivative,
            inlet_pressure_pa: pressure_derivative,
            authored_external_power_w,
            mean_stored_energy_j,
            interface_power_residual_w,
            power_residual_w: 0.0,
            solve_relative_residual: relative,
        };
        result.power_residual_w =
            result.total_dissipation_derivative_w() + terminal_power_w - authored_external_power_w;
        if ![flow_derivative, pressure_derivative]
            .iter()
            .all(|&value| complex_finite(value))
            || ![
                mean_stored_energy_j,
                authored_external_power_w,
                terminal_power_w,
                interface_power_residual_w,
                result.power_residual_w,
                result.total_dissipation_derivative_w(),
            ]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(invalid("complete DG frequency physical derivative overflows").into());
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

fn observe_terminal(
    binding: &operator::DgBinding<'_>,
    appended: &operator::AppendedDg,
    response: &DgLoadedModalHarmonicResponse,
    state: &[C64],
    derivative: &[C64],
    relative_tolerance: f64,
    gate: &CancelGate,
) -> Result<Option<DgHarmonicTerminalDerivative>, DgHarmonicCouplingError> {
    let Some(source) = &appended.terminal_source else {
        return Ok(None);
    };
    let original = response
        .terminal
        .as_ref()
        .ok_or_else(|| invalid("terminal tangent lost its original response"))?;
    let (n, ports) = (binding.terminal_dimension, binding.terminal_ports.len());
    let (start, flow_start) = (appended.terminal_offset, appended.flow_offset);
    let s = C64::new(0.0, -response.modal.omega_rad_s);
    let roundoff = 64.0 * state.len() as f64 * f64::EPSILON;
    let state_scale = state[start..flow_start]
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let derivative_scale = derivative[start..flow_start]
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let flow_scale = derivative[flow_start..appended.junction_offset]
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    for i in 0..n {
        poll(Some(gate))?;
        let inertia = s * derivative[start + i];
        let history = C64::new(0.0, 1.0) * state[start + i];
        let mut residual = inertia - history;
        let mut scale = inertia.abs() + history.abs();
        let mut coefficient_a = response.modal.omega_rad_s;
        let mut coefficient_b = 0.0;
        for j in 0..n {
            residual = residual - derivative[start + j].scale(source.a[i * n + j]);
            scale += source.a[i * n + j].abs() * derivative[start + j].abs();
            coefficient_a += source.a[i * n + j].abs();
        }
        for j in 0..ports {
            residual = residual - derivative[flow_start + j].scale(source.b[i * ports + j]);
            scale += source.b[i * ports + j].abs() * derivative[flow_start + j].abs();
            coefficient_b += source.b[i * ports + j].abs();
        }
        // Closed inputs can leave an entire retained realization row exactly
        // zero. Keep its original equation audit, allowing only the joint LU's
        // roundoff at the original source-block scales, not a relative physical
        // tolerance applied to an unrelated active row.
        let numerical_scale =
            coefficient_a * derivative_scale + coefficient_b * flow_scale + state_scale;
        if !complex_finite(residual)
            || !scale.is_finite()
            || !numerical_scale.is_finite()
            || residual.abs() > relative_tolerance * scale.max(1e-300) + roundoff * numerical_scale
        {
            return Err(
                invalid("original rational terminal frequency equation exceeds tolerance").into(),
            );
        }
    }
    let mut pressure_pa = Vec::with_capacity(ports);
    let mut mean_input_power_w = Vec::with_capacity(ports);
    for i in 0..ports {
        poll(Some(gate))?;
        let p = (0..n).fold(C64::ZERO, |sum, j| {
            sum + derivative[start + j].scale(source.c[i * n + j])
        });
        let p = (0..ports).fold(p, |sum, j| {
            sum + derivative[flow_start + j].scale(source.d[i * ports + j])
        });
        let power = power_derivative(
            original.pressure_pa[i],
            original.outward_flow_m3_s[i],
            p,
            derivative[flow_start + i],
        );
        if !complex_finite(p) || !power.is_finite() {
            return Err(invalid("original terminal pressure or power slope overflows").into());
        }
        pressure_pa.push(p);
        mean_input_power_w.push(power);
    }
    Ok(Some(DgHarmonicTerminalDerivative {
        ports: binding.terminal_ports.to_vec(),
        state: derivative[start..flow_start].to_vec(),
        outward_flow_m3_s: derivative[flow_start..appended.junction_offset].to_vec(),
        pressure_pa,
        input_power_w: mean_input_power_w.iter().sum(),
        mean_input_power_w,
    }))
}

fn observe_channel(
    view: crate::bernoulli_aperture::slit_channel::harmonic::ChannelContinuousOperator<'_>,
    response: &ChannelHarmonicResponse,
    derivative: &[C64],
    pressure_derivative: C64,
    omega: f64,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<TerminalChannelHarmonicDerivative, DgHarmonicCouplingError> {
    let n = view.dimension();
    let mass = view.cell_masses_kg();
    let area = view.cell_areas_m2();
    let conductance = view.face_conductances_kg_s();
    let mut flow = C64::ZERO;
    let (mut energy, mut heat, mut maximum_residual) = (0.0, 0.0, 0.0_f64);
    for i in 0..n {
        poll(Some(gate))?;
        let v = response.velocities_m_s[i];
        let dv = derivative[i];
        let left = if i == 0 {
            C64::ZERO
        } else {
            response.velocities_m_s[i - 1]
        };
        let d_left = if i == 0 { C64::ZERO } else { derivative[i - 1] };
        let d_right = if i + 1 == n {
            C64::ZERO
        } else {
            derivative[i + 1]
        };
        let inertia = C64::new(0.0, -omega * mass[i]) * dv;
        let history = C64::new(0.0, mass[i]) * v;
        let forcing = pressure_derivative.scale(area[i]);
        let residual = inertia
            + (dv - d_left).scale(conductance[i])
            + (dv - d_right).scale(conductance[i + 1])
            - forcing
            - history;
        let scale = inertia.abs()
            + history.abs()
            + forcing.abs()
            + conductance[i] * (dv.abs() + d_left.abs())
            + conductance[i + 1] * (dv.abs() + d_right.abs());
        if !complex_finite(residual)
            || !scale.is_finite()
            || residual.abs() > tolerance * scale.max(1e-300)
        {
            return Err(
                invalid("original channel frequency tangent violates cell momentum").into(),
            );
        }
        maximum_residual = maximum_residual.max(residual.abs() / scale.max(1e-300));
        flow = flow + dv.scale(area[i]);
        energy += 0.5 * mass[i] * (v.conj() * dv).re;
        heat += conductance[i] * ((v - left).conj() * (dv - d_left)).re;
    }
    heat += conductance[n] * (response.velocities_m_s[n - 1].conj() * derivative[n - 1]).re;
    let drive_power = pressure_derivative * response.flow_m3_s.conj();
    let flow_power = response.pressure_drop_pa * flow.conj();
    let power = (drive_power + flow_power).scale(0.5);
    let power_residual = heat - power.re;
    let reactive_residual = power.im + 2.0 * response.mean_stored_energy_j + 2.0 * omega * energy;
    let scale = 0.5 * (drive_power.abs() + flow_power.abs())
        + heat.abs()
        + 2.0 * response.mean_stored_energy_j
        + 2.0 * omega * energy.abs();
    if !complex_finite(flow)
        || !complex_finite(power)
        || ![energy, heat, scale, power_residual, reactive_residual]
            .iter()
            .all(|value| value.is_finite())
        || power_residual.abs() > tolerance * scale.max(1e-300)
        || reactive_residual.abs() > tolerance * scale.max(1e-300)
    {
        return Err(invalid(
            "original channel tangent pressure work disagrees with independent heat or storage",
        )
        .into());
    }
    Ok(TerminalChannelHarmonicDerivative {
        velocities_m_s: derivative.to_vec(),
        pressure_drop_pa: pressure_derivative,
        flow_m3_s: flow,
        mean_stored_energy_j: energy,
        viscous_dissipation_w: heat,
        input_power_w: power.re,
        power_residual_w: power_residual,
        reactive_power_residual_w: reactive_residual,
        solve_relative_residual: maximum_residual,
    })
}

#[cfg(test)]
mod tests;
