//! The original modal/Maxwell/air equations joined to a complete DG generator.
//! Retaining DG energy coordinates avoids division by a bare duct resonance.

pub use super::super::dg::ModalDgError as DgHarmonicCouplingError;
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{DgBranchLoad, DgPortLoad};
use fs_duct::dg::{DgBoundary, DgDuct, DgHarmonicResponse, PortId};
use fs_vfit::mimo::load::PreparedMimoLoad;

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
    /// Known modal, material, modal-radiation and DG storage [J]. Generic
    /// terminal storage is unavailable and is not included in this value.
    pub mean_stored_energy_j: f64,
    /// Authored power minus modal, DG and terminal input powers [W].
    pub interface_power_residual_w: f64,
    /// Known loss plus terminal exchange minus authored external power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the complete original equations.
    pub solve_relative_residual: f64,
}
impl DgLoadedModalHarmonicResponse {
    /// Original mechanical/material/radiation loss plus DG wall loss and
    /// actual matched-boundary absorption [W].
    pub fn total_dissipation_w(&self) -> f64 {
        self.modal.total_dissipation_w()
            + self.duct.wall_dissipation_w
            + self.duct.molecular_dissipation_w
            + self.duct.matched_absorption_w
            + self.branches.iter().fold(0.0, |sum, branch| {
                sum + branch.duct.wall_dissipation_w + branch.duct.molecular_dissipation_w
                    + branch.duct.matched_absorption_w
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
    /// channel of a branch belongs to its actual outer outlet. Fixed finite gap
    /// channels and driven mouths still require their own retained equations.
    ///
    /// Terminal flow is an additional unknown: p=B_out^T*x=C*z+D*u and
    /// (sI-A)z=B*u. No inverse of D or isolated terminal impedance is formed;
    /// strictly proper and singular-D loads remain valid frequency problems.
    /// The dense dimension includes every original DG and terminal state, each
    /// terminal flow and one junction flow per branch, under the same 1024 cap.
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
        if !load.terminal_channels().is_empty() || load.has_driven_mouth() {
            return Err(invalid("continuous modal/DG load response requires finite gap and driven-mouth states to be retained").into());
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
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        let view = duct.continuous_operator();
        view.admit_frequency(omega_rad_s)
            .map_err(|e| duct_error(e, gate))?;
        let n = self.mode_count();
        let terminal_ports = terminal.map_or(&[][..], |(_, ports)| ports);
        if flow_column.len() != n
            || flow_column.iter().any(|b| !b.is_finite())
            || !flow_column.iter().any(|b| *b != 0.0)
            || view.ports().len() != terminal_ports.len() + 1
            || terminal.is_some_and(|(owner, ports)| owner.ports() != ports.len())
            || terminal_ports
                .iter()
                .enumerate()
                .any(|(i, id)| *id == inlet || terminal_ports[..i].contains(id))
        {
            return Err(invalid("continuous modal/DG coupling requires a complete original volume projection and one distinct terminal channel for every non-inlet port").into());
        }
        let pressure_column = &view
            .ports()
            .iter()
            .find(|(id, _)| *id == inlet)
            .ok_or_else(|| invalid("continuous DG inlet is absent from its source basis"))?
            .1;
        let branch_views: Vec<_> = branches
            .iter()
            .map(|branch| branch.duct().continuous_operator())
            .collect();
        let mut acoustic_dimension = view.dimension();
        let mut branch_offsets = Vec::with_capacity(branches.len());
        let mut branch_ports = Vec::with_capacity(branches.len());
        let mut junction_columns = Vec::with_capacity(branches.len());
        for (index, (branch, original)) in branches.iter().zip(&branch_views).enumerate() {
            poll(Some(gate))?;
            original
                .admit_frequency(omega_rad_s)
                .map_err(|e| duct_error(e, gate))?;
            let (DgBoundary::Port(branch_inlet), DgBoundary::Port(branch_outlet)) =
                (branch.duct().inlet(), branch.duct().outlet())
            else {
                return Err(
                    invalid("continuous finite branch lost its original two flow ports").into(),
                );
            };
            if original.ports().len() != 2
                || branch_inlet == branch_outlet
                || !terminal_ports.contains(&branch.port())
                || branches[..index]
                    .iter()
                    .any(|other| other.port() == branch.port())
            {
                return Err(invalid("continuous finite branches require distinct owned junctions and complete original terminal channels").into());
            }
            let bore_column = &view
                .ports()
                .iter()
                .find(|(id, _)| *id == branch.port())
                .ok_or_else(|| invalid("continuous finite branch lost its main-bore junction"))?
                .1;
            let inlet_column = &original
                .ports()
                .iter()
                .find(|(id, _)| *id == branch_inlet)
                .ok_or_else(|| invalid("continuous finite branch lost its inlet pressure trace"))?
                .1;
            branch_offsets.push(acoustic_dimension);
            branch_ports.push((branch_inlet, branch_outlet));
            junction_columns.push((bore_column.as_slice(), inlet_column.as_slice()));
            acoustic_dimension = acoustic_dimension
                .checked_add(original.dimension())
                .ok_or_else(|| invalid("continuous finite branch dimension overflow"))?;
        }
        let terminal_columns: Vec<_> = terminal_ports
            .iter()
            .map(|id| {
                let (original, offset, local_id) =
                    match branches.iter().position(|branch| branch.port() == *id) {
                        Some(index) => (
                            &branch_views[index],
                            branch_offsets[index],
                            branch_ports[index].1,
                        ),
                        None => (&view, 0, *id),
                    };
                original
                    .ports()
                    .iter()
                    .find(|(port, _)| *port == local_id)
                    .map(|(_, column)| (offset, column.as_slice()))
                    .ok_or_else(|| {
                        invalid("continuous terminal channel names an undeclared DG port")
                    })
            })
            .collect::<Result<_, _>>()?;
        let terminal_view = terminal.map(|(owner, _)| owner.continuous_operator());
        if let Some(owner) = &terminal_view {
            owner
                .admit_frequency(omega_rad_s)
                .map_err(|error| invalid(error.0))?;
        }
        let terminal_dimension = terminal_view.as_ref().map_or(0, |owner| owner.dimension());
        let extra = acoustic_dimension
            .checked_add(terminal_dimension)
            .and_then(|value| value.checked_add(terminal_ports.len()))
            .and_then(|value| value.checked_add(branches.len()))
            .ok_or_else(|| invalid("continuous terminal harmonic dimension overflow"))?;
        let mut operator = self.harmonic_operator(omega_rad_s, force, extra, maximum_work, gate)?;
        let dimension = operator.dimension;
        let offset = operator.original_dimension;
        let terminal_offset = offset + acoustic_dimension;
        let flow_offset = terminal_offset + terminal_dimension;
        let junction_offset = flow_offset + terminal_ports.len();
        let source = view.generator_per_s();
        let terminal_source = terminal_view
            .as_ref()
            .map(|owner| {
                owner
                    .state_space(16 * dimension * dimension, &mut || gate.is_requested())
                    .map_err(|error| {
                        if gate.is_requested() {
                            ModalCouplingError::Cancelled
                        } else {
                            invalid(error.0)
                        }
                    })
            })
            .transpose()?;
        let s = C64::new(0.0, -omega_rad_s);
        for i in 0..view.dimension() {
            poll(Some(gate))?;
            for j in 0..view.dimension() {
                operator.values[(offset + i) * dimension + offset + j] =
                    C64::from_re(-source[i * view.dimension() + j]);
            }
            let diagonal = (offset + i) * dimension + offset + i;
            operator.values[diagonal] = operator.values[diagonal] + s;
            for (j, &b) in flow_column.iter().enumerate() {
                // D_modal q + b B^T x = F; (s I-L)x - s B b^T q = 0.
                let coupling = pressure_column[i] * b;
                operator.values[j * dimension + offset + i] = C64::from_re(coupling);
                operator.values[(offset + i) * dimension + j] = -s.scale(coupling);
            }
        }
        for (index, original) in branch_views.iter().enumerate() {
            let branch_offset = offset + branch_offsets[index];
            let width = original.dimension();
            let source = original.generator_per_s();
            for i in 0..width {
                poll(Some(gate))?;
                for j in 0..width {
                    operator.values[(branch_offset + i) * dimension + branch_offset + j] =
                        C64::from_re(-source[i * width + j]) + if i == j { s } else { C64::ZERO };
                }
            }
            // Main outward q and branch inward q share one multiplier. Its
            // equation enforces the ORIGINAL pressure traces, not an impedance
            // condensation, including when either isolated duct is singular.
            let junction = junction_offset + index;
            for (i, &b) in junction_columns[index].0.iter().enumerate() {
                operator.values[(offset + i) * dimension + junction] = C64::from_re(b);
                operator.values[junction * dimension + offset + i] = C64::from_re(b);
            }
            for (i, &b) in junction_columns[index].1.iter().enumerate() {
                operator.values[(branch_offset + i) * dimension + junction] = C64::from_re(-b);
                operator.values[junction * dimension + branch_offset + i] = C64::from_re(-b);
            }
        }
        if let Some(load) = &terminal_source {
            let ports = terminal_ports.len();
            for (j, &(relative_offset, column)) in terminal_columns.iter().enumerate() {
                poll(Some(gate))?;
                let terminal_duct_offset = offset + relative_offset;
                for (i, &b) in column.iter().enumerate() {
                    operator.values[(terminal_duct_offset + i) * dimension + flow_offset + j] =
                        C64::from_re(b);
                    operator.values[(flow_offset + j) * dimension + terminal_duct_offset + i] =
                        C64::from_re(b);
                }
            }
            for i in 0..terminal_dimension {
                poll(Some(gate))?;
                for j in 0..terminal_dimension {
                    operator.values[(terminal_offset + i) * dimension + terminal_offset + j] =
                        C64::from_re(-load.a[i * terminal_dimension + j])
                            + if i == j { s } else { C64::ZERO };
                }
                for j in 0..ports {
                    operator.values[(terminal_offset + i) * dimension + flow_offset + j] =
                        C64::from_re(-load.b[i * ports + j]);
                    operator.values[(flow_offset + j) * dimension + terminal_offset + i] =
                        C64::from_re(-load.c[j * terminal_dimension + i]);
                }
            }
            for i in 0..ports {
                for j in 0..ports {
                    operator.values[(flow_offset + i) * dimension + flow_offset + j] =
                        C64::from_re(-load.d[i * ports + j]);
                }
            }
        }
        let (state, relative) = solve_harmonic(
            &operator.values,
            dimension,
            force,
            self.config.solve_relative_tolerance,
            gate,
        )?;
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
        for (index, (branch, original)) in branches.iter().zip(&branch_views).enumerate() {
            poll(Some(gate))?;
            let local_offset = offset + branch_offsets[index];
            let terminal_index = terminal_ports
                .iter()
                .position(|id| *id == branch.port())
                .ok_or_else(|| {
                    invalid("continuous branch observation lost its terminal channel")
                })?;
            let junction_flow = state[junction_offset + index];
            let inputs = [
                (branch_ports[index].0, junction_flow),
                (branch_ports[index].1, -state[flow_offset + terminal_index]),
            ];
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
        );
        let interface_power_residual_w = authored_external_power_w
            - modal.external_power_w
            - duct.input_power_w
            - branch_responses
                .iter()
                .map(|branch| branch.duct.input_power_w)
                .sum::<f64>()
            - terminal_power_w;
        let mut result = DgLoadedModalHarmonicResponse {
            modal,
            duct,
            branches: branch_responses,
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
