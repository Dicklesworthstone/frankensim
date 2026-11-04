//! The original modal/Maxwell/air equations joined to a complete DG generator.
//! Retaining DG energy coordinates avoids division by a bare duct resonance.

pub use super::super::dg::ModalDgError as DgHarmonicCouplingError;
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{DgBranchLoad, DgPortLoad, DgTerminalChannelKind, DgTerminalChannelLoad};
use crate::bernoulli_aperture::slit_channel::harmonic::ChannelHarmonicResponse;
use fs_duct::dg::{DgBoundary, DgDuct, DgHarmonicResponse, PortId};
use fs_vfit::mimo::load::PreparedMimoLoad;

mod operator;
mod frequency;
pub use frequency::{DgHarmonicBranchDerivative, DgHarmonicChannelDerivative,
    DgHarmonicTerminalDerivative, DgLoadedModalHarmonicFrequencyResponse,
    TerminalChannelHarmonicDerivative};
/// Free modes of the full original mechanical, material and acoustic pencil.
pub mod poles;
/// Driven physical peaks and widths with every attached acoustic owner retained.
pub mod resonance;

/// Original terminal real-state phasors and physical outward port exchange.
/// The admitted source's passivity reports do not supply a generic storage
/// matrix, so no terminal stored energy or heat is invented here.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicTerminalResponse {
    /// Original explicit terminal channel order, including mutual entries.
    pub ports: Vec<PortId>,
    /// Peak phasors of the original REAL state-space coordinates. A conjugate
    /// pair uses [Re(z) for all inputs, -Im(z) for all inputs].
    pub state: Vec<C64>,
    /// Outward flow from each original duct port into its terminal [m3/s].
    pub outward_flow_m3_s: Vec<C64>,
    /// Original C*z+D*u pressure in each terminal channel [Pa].
    pub pressure_pa: Vec<C64>,
    /// Signed mean work entering each original terminal channel [W].
    pub mean_input_power_w: Vec<f64>,
    /// Sum of the signed channel work rates [W].
    pub input_power_w: f64,
}

/// One original finite branch and both sides of its zero-storage junction.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicBranchResponse {
    /// Original main-bore side-port ID. The same external terminal channel
    /// observes this branch's OUTLET, not its internal inlet junction.
    pub port: PortId,
    /// Positive from the main bore into the finite branch [m3/s].
    pub junction_flow_m3_s: C64,
    /// Original main-bore pressure at this internal junction [Pa].
    pub bore_pressure_pa: C64,
    /// Independently observed original branch inlet pressure [Pa].
    pub inlet_pressure_pa: C64,
    /// Signed cycle-mean power OUT of the main bore [W].
    pub bore_power_w: f64,
    /// Signed cycle-mean power INTO the branch [W].
    pub branch_power_w: f64,
    /// Bore power minus branch power [W]; a junction creates no storage or loss.
    pub junction_power_residual_w: f64,
    /// Complete original branch coordinates, fields, local ports and losses.
    pub duct: DgHarmonicResponse,
}

/// Original finite terminal passage or rigid closure in the joint AC solve.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicChannelResponse {
    pub port: PortId,
    pub duct_pressure_pa: C64,
    pub exterior_pressure_pa: C64,
    pub outward_flow_m3_s: C64,
    /// Original pressure-drop multiplier, or rigid closure reaction [Pa].
    pub pressure_drop_pa: C64,
    /// Observed duct minus exterior pressure minus multiplier [Pa].
    pub pressure_drop_residual_pa: C64,
    /// Original field volume flow minus the shared terminal flow [m3/s].
    pub flow_residual_m3_s: C64,
    /// Observed pressure-plane difference times shared flow, cycle mean [W].
    pub input_power_w: f64,
    /// Shared pressure work minus the original cell field's pressure work [W].
    pub interface_power_residual_w: f64,
    /// Every original cell and its physical storage/heat; None for a closure.
    pub channel: Option<ChannelHarmonicResponse>,
}

/// Continuous AC response of the original connected owners, in exp(-i omega t).
#[derive(Clone, Debug, PartialEq)]
pub struct DgLoadedModalHarmonicResponse {
    /// Original modal response. Its input power includes the opposing duct
    /// reaction; `authored_external_power_w` is the external power of the joint
    /// system. Its residual describes the complete coupled equations.
    pub modal: ModalHarmonicResponse,
    /// Actual DG energy coordinates, physical fields, port pressures and losses.
    pub duct: DgHarmonicResponse,
    /// Every original finite side branch, in its explicit attachment order.
    pub branches: Vec<DgHarmonicBranchResponse>,
    /// Original finite terminal fields and rigid closures, in source order.
    pub channels: Vec<DgHarmonicChannelResponse>,
    /// Complete original rational radiation/terminal response, when attached.
    pub terminal: Option<DgHarmonicTerminalResponse>,
    /// Signed mean exchange into the terminal, separate from known losses [W].
    pub terminal_power_w: f64,
    /// Physical flow INTO the original duct inlet [m3/s].
    pub inlet_flow_m3_s: C64,
    /// Physical inlet pressure opposing that flow [Pa].
    pub inlet_pressure_pa: C64,
    /// Mean power of the authored generalized forces only [W].
    pub authored_external_power_w: f64,
    /// Known modal, material, modal-radiation, DG and finite-channel storage [J]. Generic
    /// terminal storage is unavailable and is not included in this value.
    pub mean_stored_energy_j: f64,
    /// Authored power minus modal, DG, finite-channel and terminal input powers [W].
    pub interface_power_residual_w: f64,
    /// Known loss plus terminal exchange minus authored external power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the complete original equations.
    pub solve_relative_residual: f64,
}
impl DgLoadedModalHarmonicResponse {
    /// Original mechanical/material/radiation loss, DG wall loss, actual
    /// matched-boundary absorption and finite-channel no-slip wall heat [W].
    pub fn total_dissipation_w(&self) -> f64 {
        self.modal.total_dissipation_w()
            + self.duct.wall_dissipation_w
            + self.duct.molecular_dissipation_w
            + self.duct.matched_absorption_w
            + self.branches.iter().fold(0.0, |sum, branch| {
                sum + branch.duct.wall_dissipation_w + branch.duct.molecular_dissipation_w
                    + branch.duct.matched_absorption_w
            })
            + self.channels.iter().fold(0.0, |sum, channel| {
                sum + channel.channel.as_ref().map_or(0.0, |field| field.viscous_dissipation_w)
            })
    }
}

fn duct_error(error: fs_duct::DuctError, gate: &CancelGate) -> DgHarmonicCouplingError {
    if gate.is_requested() {
        ModalCouplingError::Cancelled.into()
    } else {
        DgHarmonicCouplingError::Acoustic(AcousticRealizeError::Duct(error))
    }
}

impl CoupledModalSystem {
    /// Continuous forced response with one complete source-bound DG duct.
    ///
    /// `flow_column` maps ORIGINAL modal velocities to inlet volume flow:
    /// Q=b^T v. The opposite generalized force is -b p. The duct retains all
    /// its original energy coordinates, including wall memory and its actual
    /// closed/open/matched distal boundary. Exactly one declared DG flow port
    /// must exist and equal `inlet`; omitted external ports are never closed
    /// implicitly. Branches and rational terminations require their own retained
    /// states and are not supplied by a bare DgDuct.
    ///
    /// The existing complex LU solves the complete matrix, without a standalone
    /// impedance inversion, fitted poles or a coupling timestep. Constitutive
    /// band and spatial-resolution obligations stay with their original owners.
    /// No accepted history changes. The dense allowance is
    /// d^3+d^2*(connections+material_arms+5), with
    /// d=original_mechanical_modes+air_poles+DG_states <=1024.
    #[allow(clippy::too_many_arguments)]
    pub fn harmonic_response_with_duct(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        self.harmonic_response_bound(
            omega_rad_s,
            force,
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

    /// Retain the original DG bore, finite branches and COMPLETE terminal matrix.
    /// Every explicit terminal channel, pole state and mutual coefficient is
    /// retained in original order. The terminal's original validity band applies.
    /// Each finite branch keeps its entire source operator and its own material
    /// band. Junction flows and pressure continuity are solved simultaneously;
    /// no branch impedance is inverted at a bare resonant pole. The terminal
    /// channel of a branch belongs to its actual outer outlet. Each finite gap
    /// retains every original transverse/radial momentum cell and one pressure
    /// drop; rigid closures enforce zero through-flow without deleting exterior
    /// states. This is a zero-mean perturbation of the original fixed geometry
    /// and source gas; the declared resolved-speed Mach limit still applies.
    /// Driven mouths require their own retained source equations.
    ///
    /// Terminal flow is an additional unknown: p=B_out^T*x=C*z+D*u and
    /// (sI-A)z=B*u. No inverse of D or isolated terminal impedance is formed;
    /// strictly proper and singular-D loads remain valid frequency problems.
    /// The dense dimension includes every original DG and terminal state, each
    /// terminal flow, one junction flow per branch, and all finite-gap cells
    /// plus their pressure drops, under the same 1024 cap and work allowance.
    /// Original channel observations additionally reserve 256*f+64*d*c work,
    /// where f counts finite cells plus their drops and c counts channel inputs.
    /// Both reservations precede the joint dense and returned-field allocations.
    /// Accepted states and source/passivity evidence remain unchanged. Terminal
    /// work is reported separately; an admitted fit supplies no generic energy
    /// certificate for its internal coordinates.
    pub fn harmonic_response_with_load(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        load: &DgPortLoad,
        flow_column: &[f64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        if load.has_driven_mouth() {
            return Err(invalid("continuous modal/DG load response requires driven-mouth source equations to be retained").into());
        }
        let terminal = load.terminal().map(|owner| (owner, load.terminal_ports()));
        self.harmonic_response_bound(
            omega_rad_s,
            force,
            load.duct(),
            load.inlet(),
            flow_column,
            terminal,
            load.branches(),
            load.terminal_channels(),
            maximum_work,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn harmonic_response_bound(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        terminal: Option<(&PreparedMimoLoad, &[PortId])>,
        branches: &[DgBranchLoad],
        channels: &[DgTerminalChannelLoad],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        let view = duct.continuous_operator();
        let n = self.mode_count();
        let branch_views: Vec<_> = branches.iter()
            .map(|branch| branch.duct().continuous_operator()).collect();
        let binding = operator::DgBinding::new(n, &view, inlet, flow_column,
            terminal, branches, &branch_views, channels, gate)?;
        binding.admit_frequency(omega_rad_s, gate)?;
        let channel_dimension = binding.channel_dimension;
        let extra = binding.additional_dimension()?;
        let full_dimension = n.checked_add(self.radiation.as_ref().map_or(0, |r| r.model.poles.len()))
            .and_then(|value| value.checked_add(extra))
            .ok_or_else(|| invalid("continuous channel joint dimension overflow"))?;
        let channel_work = channel_dimension.checked_mul(256)
            .and_then(|value| channels.len().checked_mul(full_dimension)
                .and_then(|value| value.checked_mul(64))
                .and_then(|interface| value.checked_add(interface)))
            .ok_or_else(|| invalid("continuous channel observation work overflow"))?;
        let operator_work = maximum_work.checked_sub(channel_work)
            .ok_or_else(|| invalid("continuous channel observations exceed the declared joint work allowance"))?;
        let mut operator = self.harmonic_operator(omega_rad_s, force, extra, operator_work, gate)?;
        let dimension = operator.dimension;
        let offset = operator.original_dimension;
        let s = C64::new(0.0, -omega_rad_s);
        let appended = binding.append(offset, flow_column,
                16 * dimension * dimension,
                &mut operator::AcEntries { values: &mut operator.values, dimension, s }, gate)?;
        let (state, relative) = solve_harmonic(
            &operator.values,
            dimension,
            force,
            self.config.solve_relative_tolerance,
            gate,
        )?;
        self.observe_harmonic_dg_state(omega_rad_s, force, inlet, flow_column, branches,
            &binding, &appended, &state, relative, maximum_work, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_harmonic_dg_state(
        &self, omega_rad_s: f64, force: &[C64], inlet: PortId, flow_column: &[f64],
        branches: &[DgBranchLoad], binding: &operator::DgBinding<'_>,
        appended: &operator::AppendedDg, state: &[C64], relative: f64,
        maximum_work: usize, gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        let n = self.mode_count();
        let dimension = state.len();
        let offset = dimension - binding.additional_dimension()?;
        let s = C64::new(0.0, -omega_rad_s);
        let view = binding.view;
        let branch_views = binding.branch_views;
        let terminal_ports = binding.terminal_ports;
        let terminal_dimension = binding.terminal_dimension;
        let branch_offsets = &binding.branch_offsets;
        let branch_ports = &binding.branch_ports;
        let junction_columns = &binding.junction_columns;
        let terminal_columns = &binding.terminal_columns;
        let channels = binding.channels;
        let terminal_source = &appended.terminal_source;
        let terminal_offset = appended.terminal_offset;
        let flow_offset = appended.flow_offset;
        let junction_offset = appended.junction_offset;
        let channel_offsets = &appended.channel_offsets;
        let flow = s * projection(flow_column, &state[..n]);
        let mut inflows = Vec::with_capacity(terminal_ports.len() + 1);
        inflows.push((inlet, flow));
        for (j, &id) in terminal_ports.iter().enumerate() {
            let outward = match branches.iter().position(|branch| branch.port() == id) {
                Some(index) => state[junction_offset + index],
                None => state[flow_offset + j],
            };
            inflows.push((id, -outward));
        }
        for (index, branch) in branches.iter().enumerate() {
            if branch_ports[index].1.is_none() {
                inflows.push((branch.port(), -state[junction_offset + index]));
            }
        }
        let terminal = terminal_source
            .as_ref()
            .map(|load| {
                let ports = terminal_ports.len();
                let mut pressure_pa = vec![C64::ZERO; ports];
                let mut mean_input_power_w = Vec::with_capacity(ports);
                for i in 0..ports {
                    poll(Some(gate))?;
                    let p = (0..terminal_dimension).fold(C64::ZERO, |sum, j| {
                        sum + state[terminal_offset + j].scale(load.c[i * terminal_dimension + j])
                    });
                    pressure_pa[i] = (0..ports).fold(p, |sum, j| {
                        sum + state[flow_offset + j].scale(load.d[i * ports + j])
                    });
                    let power = 0.5 * (pressure_pa[i] * state[flow_offset + i].conj()).re;
                    if !complex_finite(pressure_pa[i]) || !power.is_finite() {
                        return Err(invalid("continuous terminal pressure or work overflows"));
                    }
                    mean_input_power_w.push(power);
                }
                Ok(DgHarmonicTerminalResponse {
                    ports: terminal_ports.to_vec(),
                    state: state[terminal_offset..flow_offset].to_vec(),
                    outward_flow_m3_s: state[flow_offset..junction_offset].to_vec(),
                    pressure_pa,
                    input_power_w: mean_input_power_w.iter().sum(),
                    mean_input_power_w,
                })
            })
            .transpose()?;
        let terminal_power_w = terminal.as_ref().map_or(0.0, |value| value.input_power_w);
        let mut cancelled = || gate.is_requested();
        let duct = view
            .observe_harmonic(
                omega_rad_s,
                &state[offset..offset + view.dimension()],
                &inflows,
                maximum_work,
                self.config.solve_relative_tolerance,
                &mut cancelled,
            )
            .map_err(|e| duct_error(e, gate))?;
        let mut branch_responses = Vec::with_capacity(branches.len());
        for (index, (branch, original)) in branches.iter().zip(branch_views).enumerate() {
            poll(Some(gate))?;
            let local_offset = offset + branch_offsets[index];
            let junction_flow = state[junction_offset + index];
            let mut inputs = vec![(branch_ports[index].0, junction_flow)];
            if let Some(outlet) = branch_ports[index].1 {
                let terminal_index = terminal_ports.iter()
                    .position(|id| *id == branch.port())
                    .ok_or_else(|| invalid("continuous branch observation lost its terminal channel"))?;
                inputs.push((outlet, -state[flow_offset + terminal_index]));
            }
            let response = original
                .observe_harmonic(
                    omega_rad_s,
                    &state[local_offset..local_offset + original.dimension()],
                    &inputs,
                    maximum_work,
                    self.config.solve_relative_tolerance,
                    &mut cancelled,
                )
                .map_err(|e| duct_error(e, gate))?;
            let bore_pressure_pa = duct
                .ports
                .iter()
                .find(|port| port.id == branch.port())
                .ok_or_else(|| invalid("continuous bore observation omitted a branch junction"))?
                .pressure_pa;
            let inlet_pressure_pa = response
                .ports
                .iter()
                .find(|port| port.id == branch_ports[index].0)
                .ok_or_else(|| invalid("continuous branch observation omitted its inlet"))?
                .pressure_pa;
            let bore_power_w = 0.5 * (bore_pressure_pa * junction_flow.conj()).re;
            let branch_power_w = 0.5 * (inlet_pressure_pa * junction_flow.conj()).re;
            let junction_power_residual_w = bore_power_w - branch_power_w;
            // A global backward norm can underweight this pressure row when
            // other operators carry much larger units. Check both original
            // pressure traces with their uncancelled physical terms, then the
            // independently observed exchange; neither gate borrows that norm.
            let pressure_scale = junction_columns[index]
                .0
                .iter()
                .enumerate()
                .map(|(i, &b)| b.abs() * state[offset + i].abs())
                .sum::<f64>()
                + junction_columns[index]
                    .1
                    .iter()
                    .enumerate()
                    .map(|(i, &b)| b.abs() * state[local_offset + i].abs())
                    .sum::<f64>();
            let relative = self.config.solve_relative_tolerance
                + 64.0 * f64::EPSILON * (view.dimension() + original.dimension()) as f64;
            let power_scale = 0.5 * pressure_scale * junction_flow.abs()
                + bore_power_w.abs()
                + branch_power_w.abs();
            if ![
                bore_power_w,
                branch_power_w,
                junction_power_residual_w,
                pressure_scale,
                power_scale,
            ]
            .iter()
            .all(|power| power.is_finite())
            {
                return Err(invalid("continuous branch junction power overflows").into());
            }
            if (bore_pressure_pa - inlet_pressure_pa).abs() > relative * pressure_scale.max(1e-300)
                || junction_power_residual_w.abs() > relative * power_scale.max(1e-300)
            {
                return Err(invalid("continuous original bore and finite branch disagree on junction pressure or power").into());
            }
            branch_responses.push(DgHarmonicBranchResponse {
                port: branch.port(),
                junction_flow_m3_s: junction_flow,
                bore_pressure_pa,
                inlet_pressure_pa,
                bore_power_w,
                branch_power_w,
                junction_power_residual_w,
                duct: response,
            });
        }
        let mut channel_responses = Vec::with_capacity(channels.len());
        for (index, owner) in channels.iter().enumerate() {
            poll(Some(gate))?;
            let port = owner.input().port;
            let port_index = terminal_ports.iter().position(|id| *id == port)
                .ok_or_else(|| invalid("continuous channel observation lost its terminal row"))?;
            let exterior = terminal.as_ref()
                .ok_or_else(|| invalid("continuous channel lost its complete exterior response"))?;
            let exterior_source = terminal_source.as_ref()
                .ok_or_else(|| invalid("continuous channel lost its original exterior equations"))?;
            let (observed, observed_port) = match branches.iter().position(|branch| branch.port() == port) {
                Some(branch) => (&branch_responses[branch].duct, branch_ports[branch].1
                    .ok_or_else(|| invalid("continuous finite gap cannot follow a sealed branch"))?),
                None => (&duct, port),
            };
            let duct_pressure = observed.ports.iter().find(|value| value.id == observed_port)
                .ok_or_else(|| invalid("continuous channel lost its physical pressure-plane observation"))?
                .pressure_pa;
            let exterior_pressure = exterior.pressure_pa[port_index];
            let drop = duct_pressure - exterior_pressure;
            let outward = exterior.outward_flow_m3_s[port_index];
            let input_power = 0.5 * (drop * outward.conj()).re;
            let (duct_offset, column) = terminal_columns[port_index];
            let pressure_scale = column.iter().enumerate().map(|(i, b)|
                b.abs() * state[offset + duct_offset + i].abs()).sum::<f64>()
                + (0..terminal_dimension).map(|i|
                    exterior_source.c[port_index * terminal_dimension + i].abs()
                        * state[terminal_offset + i].abs()).sum::<f64>()
                + (0..terminal_ports.len()).map(|i|
                    exterior_source.d[port_index * terminal_ports.len() + i].abs()
                        * state[flow_offset + i].abs()).sum::<f64>();
            let (constitutive_drop, flow_residual, field, flow_scale, pressure_residual) =
                if let Some(original) = owner.channel() {
                    let view = original.continuous_operator()?;
                    let start = channel_offsets[index];
                    let velocities = &state[start..start + view.dimension()];
                    let constitutive_drop = state[start + view.dimension()];
                    let field = view.observe_harmonic(
                        omega_rad_s, velocities, drop, maximum_work,
                        self.config.solve_relative_tolerance, &mut || gate.is_requested(),
                    )?;
                    let maximum_mach = match owner.input().kind {
                        DgTerminalChannelKind::Slit { maximum_mach, .. }
                        | DgTerminalChannelKind::Circular { maximum_mach, .. } => maximum_mach,
                        DgTerminalChannelKind::Closed => return Err(invalid(
                            "continuous finite channel lost its physical geometry").into()),
                    };
                    if velocities.iter().any(|v| v.abs() > maximum_mach * owner.gas().sound_speed) {
                        return Err(invalid("continuous terminal channel exceeds its declared local resolved-speed Mach limit").into());
                    }
                    let flow_scale = outward.abs() + view.cell_areas_m2().iter()
                        .zip(velocities).map(|(area, v)| area * v.abs()).sum::<f64>();
                    (constitutive_drop, field.flow_m3_s - outward, Some(field),
                        flow_scale, drop - constitutive_drop)
                } else {
                    // Zero-through-flow closure retains the exterior reaction
                    // and mutual states, with no fabricated fluid storage.
                    (drop, -outward, None, flow.abs() + exterior.outward_flow_m3_s.iter()
                        .map(|value| value.abs()).sum::<f64>(), C64::ZERO)
                };
            let pressure_scale = pressure_scale + constitutive_drop.abs();
            let field_power = field.as_ref().map_or(0.0, |value| value.input_power_w);
            let work_residual = input_power - field_power;
            let work_scale = input_power.abs() + field_power.abs()
                + 0.5 * pressure_scale * flow_scale;
            let tolerance = self.config.solve_relative_tolerance
                + 64.0 * f64::EPSILON * dimension as f64;
            if ![drop, constitutive_drop, outward, flow_residual, pressure_residual]
                .iter().all(|&v| complex_finite(v))
                || ![input_power, flow_scale, pressure_scale, work_scale, work_residual]
                    .iter().all(|value| value.is_finite())
                || pressure_residual.abs() > tolerance * pressure_scale.max(1e-300)
                || flow_residual.abs() > tolerance * flow_scale.max(1e-300)
                || work_residual.abs() > tolerance * work_scale.max(1e-300)
            {
                return Err(invalid("continuous original terminal field disagrees with shared pressure, flow or work").into());
            }
            channel_responses.push(DgHarmonicChannelResponse {
                port, duct_pressure_pa: duct_pressure, exterior_pressure_pa: exterior_pressure,
                outward_flow_m3_s: outward, pressure_drop_pa: constitutive_drop,
                pressure_drop_residual_pa: pressure_residual, flow_residual_m3_s: flow_residual,
                input_power_w: input_power, interface_power_residual_w: work_residual, channel: field,
            });
        }
        let channel_power_w = channel_responses.iter().map(|channel| channel.input_power_w).sum::<f64>();
        let channel_storage_j = channel_responses.iter().map(|channel|
            channel.channel.as_ref().map_or(0.0, |field| field.mean_stored_energy_j)).sum::<f64>();
        let pressure = duct
            .ports
            .iter()
            .find(|p| p.id == inlet)
            .ok_or_else(|| invalid("continuous DG observation omitted its original inlet"))?
            .pressure_pa;
        let applied: Vec<_> = force
            .iter()
            .zip(flow_column)
            .map(|(&f, &b)| f - pressure.scale(b))
            .collect();
        if applied.iter().any(|&f| !complex_finite(f)) {
            return Err(invalid("continuous DG reaction force overflows").into());
        }
        let modal =
            self.observe_harmonic_state(omega_rad_s, &applied, &state[..offset], relative, gate)?;
        let authored_external_power_w = force
            .iter()
            .zip(&modal.velocity_m_sqrt_kg_per_s)
            .map(|(&f, &v)| 0.5 * (f * v.conj()).re)
            .sum::<f64>();
        let mean_stored_energy_j = branch_responses.iter().fold(
            modal.mean_stored_energy_j + duct.mean_stored_energy_j,
            |sum, branch| sum + branch.duct.mean_stored_energy_j,
        ) + channel_storage_j;
        let interface_power_residual_w = authored_external_power_w
            - modal.external_power_w
            - duct.input_power_w
            - branch_responses
                .iter()
                .map(|branch| branch.duct.input_power_w)
                .sum::<f64>()
            - terminal_power_w
            - channel_power_w;
        let mut result = DgLoadedModalHarmonicResponse {
            modal,
            duct,
            branches: branch_responses,
            channels: channel_responses,
            terminal,
            terminal_power_w,
            inlet_flow_m3_s: flow,
            inlet_pressure_pa: pressure,
            authored_external_power_w,
            mean_stored_energy_j,
            interface_power_residual_w,
            power_residual_w: 0.0,
            solve_relative_residual: relative,
        };
        result.power_residual_w =
            result.total_dissipation_w() + terminal_power_w - authored_external_power_w;
        if [
            result.authored_external_power_w,
            result.mean_stored_energy_j,
            result.interface_power_residual_w,
            result.power_residual_w,
            result.total_dissipation_w(),
            result.terminal_power_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("continuous modal/DG physical power or storage overflows").into());
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod terminal_tests;

#[cfg(test)]
mod mimo_tests;

#[cfg(test)]
mod branch_tests;

