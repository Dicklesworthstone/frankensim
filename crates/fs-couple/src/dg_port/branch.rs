//! Finite side ducts on the original bore's work-conjugate interface ports.
//!
//! Every branch keeps its original DG geometry, gas, wall memory and exact
//! held-flow propagator. Only zero-storage junction flows are eliminated.
//! This is a point junction plus a resolved finite duct, not a three-dimensional
//! entrance flow model. Junction-clock and spatial refinement remain necessary.

use super::*;
use fs_duct::dg::PortOutput;
use fs_material::gas::GasState;

mod coupling;
pub(super) use coupling::BranchCoupling;

/// One original finite duct attached to an existing bore side port.
#[derive(Debug, Clone)]
pub struct DgBranch {
    /// Stable external channel ID. After attachment it denotes the branch's
    /// OUTER outlet, while the same original bore ID is its internal junction.
    pub port: PortId,
    /// One inlet flow port and either an outlet flow port or an exact rigid
    /// closure. A sealed branch retains all acoustic and wall-memory state but
    /// has no exterior radiation channel. Local PortIds belong to this duct.
    /// Side ports on a branch are refused by this bounded star topology.
    pub duct: DgDuct,
}

/// Original retained branch state and propagation owner; no mutable escape.
#[derive(Debug, Clone)]
pub struct DgBranchLoad {
    port: PortId,
    duct: DgDuct,
    propagator: DgPropagator,
    pub(super) state: DgState,
    inlet: PortId,
    outlet: Option<PortId>,
}
impl DgBranchLoad {
    #[must_use]
    pub const fn port(&self) -> PortId { self.port }
    #[must_use]
    pub const fn duct(&self) -> &DgDuct { &self.duct }
    #[must_use]
    pub const fn propagator(&self) -> &DgPropagator { &self.propagator }
    #[must_use]
    pub const fn state(&self) -> &DgState { &self.state }

    pub(super) fn prepare_gas(
        &self,
        field: &fs_material::gas_field::AxialGasField,
        refresh: DgGasRefresh,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, fs_duct::dg::DgGasMaterialUpdate), AcousticRealizeError> {
        let (duct, state, propagator, report) = prepare_dg_gas(
            &self.duct, &self.state, &self.propagator, field, refresh,
            maximum_absolute_parameter_work_j, cancelled,
        )?;
        Ok((Self {
            port: self.port, duct, propagator, state,
            inlet: self.inlet, outlet: self.outlet,
        }, report))
    }

    pub(super) fn prepare(
        source: DgBranch,
        dt: f64,
        tolerance: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(inlet) = source.duct.inlet() else {
            return Err(invalid(
                "DG finite branch requires its original inlet flow port",
            ));
        };
        let outlet = match source.duct.outlet() {
            DgBoundary::Port(port) => Some(port),
            DgBoundary::Closed => None,
            _ => {
                return Err(invalid(
                    "DG finite branch requires an outlet flow port or exact closure",
                ));
            }
        };
        let propagator = source.duct.propagator(dt, tolerance)
            .map_err(AcousticRealizeError::Duct)?;
        let state = source.duct.zero_state();
        Ok(Self { port: source.port, duct: source.duct, propagator, state, inlet, outlet })
    }
}

/// Both independently integrated sides of one accepted internal junction.
#[derive(Debug, Clone, PartialEq)]
pub struct DgBranchFrame {
    pub port: PortId,
    /// Positive from the main bore INTO the finite branch [m³/s]. This is not
    /// the exterior flow: compressibility and memory separate the two ends.
    pub junction_flow_m3_s: f64,
    pub bore_port: PortOutput,
    pub branch_inlet: PortOutput,
    /// Signed work OUT of the bore at its junction [J].
    pub bore_work_j: f64,
    /// Signed work INTO the branch at its inlet [J].
    pub branch_work_j: f64,
    /// Bore work minus branch work [J]; checked before any history is accepted.
    pub junction_work_residual_j: f64,
    /// Unmodified original branch step, including its local IDs and ledger.
    pub duct: DgStepOutput,
}

impl DgPortLoad {
    /// Attach sealed finite branches to a bore with no radiating flow ports.
    /// All branch volumes and original propagation states remain in the source
    /// load, including their same-step junction feedback and signed work.
    pub fn new_with_closed_branches(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        branches: Vec<DgBranch>,
    ) -> Result<Self, AcousticRealizeError> {
        Self::prepare(
            duct,
            time_step_s,
            exponential_tolerance,
            maximum_dofs,
            Vec::new(),
            None,
            branches,
        )
    }

    /// Attach finite original DG ducts before the complete exterior MIMO load.
    ///
    /// The original bore side-port IDs remain the external load channel IDs;
    /// their pressures and flows now belong to the branch OUTLETS. Every bore
    /// unsealed port except the source inlet is explicitly owned by the load.
    /// Sealed branch junctions remain internal and consume no exterior row.
    /// Unattached side ports remain flush point apertures.
    ///
    /// All ducts advance on one clock. The complete acoustic and wall-state
    /// sum must fit maximum_dofs before any dense exponential is prepared.
    /// Junction mean pressures agree and opposite flow work cancels before
    /// publication. This held-flow interconnection is not the joint continuous
    /// exponential and does not add a finite-volume junction or entrance loss.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_terminal_branches(
        duct: DgDuct, time_step_s: f64, exponential_tolerance: f64,
        maximum_dofs: usize, ports: Vec<PortId>, terminal: PreparedMimoLoad,
        branches: Vec<DgBranch>,
    ) -> Result<Self, AcousticRealizeError> {
        Self::prepare(duct, time_step_s, exponential_tolerance, maximum_dofs,
            ports, Some(terminal), branches)
    }

    /// All original finite branches, in their explicit attachment order.
    #[must_use]
    pub fn branches(&self) -> &[DgBranchLoad] { &self.branches }

    /// Accepted internal work exchanges; empty until the first accepted step.
    #[must_use]
    pub fn last_branch_frames(&self) -> &[DgBranchFrame] { &self.last_branch_frames }

    /// Complete bore/branch acoustic, wall-memory and finite gap kinetic storage [J].
    /// Original rational terminal storage is not fabricated by this sum.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        let ducts = self.branches.iter().fold(self.state.stored_energy_j(),
            |sum, branch| sum + branch.state.stored_energy_j());
        self.terminal_channels.iter().fold(ducts, |sum, channel| sum + channel.stored_energy_j())
    }

    /// Original inferred duct loss plus gap wall heat and BE numerical loss [J].
    /// Gap frames/cumulative getters retain those two distinct contributions.
    #[must_use]
    pub fn cumulative_dissipated_j(&self) -> f64 {
        let ducts = self.branches.iter().fold(self.state.cumulative_dissipated_j(),
            |sum, branch| sum + branch.state.cumulative_dissipated_j());
        self.terminal_channels.iter().fold(ducts, |sum, channel| sum + channel.cumulative_dissipated_j())
    }

    /// Tightest declared viscothermal upper frequency across every actual duct.
    #[must_use]
    pub fn minimum_wall_band_upper_rad_s(&self) -> Option<f64> {
        std::iter::once(&self.duct).chain(self.branches.iter().map(|b| &b.duct))
            .filter_map(|duct| duct.viscothermal_band_rad_s().map(|band| band[1]))
            .reduce(f64::min)
    }

    /// Material gas at the actual external load channel, after any finite branch.
    pub fn terminal_gas(&self, port: PortId) -> Result<GasState, AcousticRealizeError> {
        if !self.terminal_ports().contains(&port) {
            return Err(invalid("undeclared DG terminal material port"));
        }
        if let Some(branch) = self.branches.iter().find(|branch| branch.port == port) {
            return Ok(*branch.duct.outlet_gas());
        }
        if self.duct.outlet() == DgBoundary::Port(port) {
            return Ok(*self.duct.outlet_gas());
        }
        self.duct.side_port_source(port).map(|source| source.gas)
            .map_err(|_| invalid("DG terminal lost its original local side-port gas"))
    }

    /// Actual outer duct area, when geometry specifies one. A flush point side
    /// port has no independently defined aperture area and returns None.
    pub fn terminal_area_m2(&self, port: PortId) -> Result<Option<f64>, AcousticRealizeError> {
        if !self.terminal_ports().contains(&port) {
            return Err(invalid("undeclared DG terminal geometry port"));
        }
        if let Some(branch) = self.branches.iter().find(|branch| branch.port == port) {
            return Ok(Some(branch.duct.outlet_area_m2()));
        }
        Ok((self.duct.outlet() == DgBoundary::Port(port)).then(|| self.duct.outlet_area_m2()))
    }

    pub(super) fn multiport_response(&self) -> Result<DgMultiportResponse, AcousticRealizeError> {
        let main = self.propagator.multiport_response(&self.state)
            .map_err(AcousticRealizeError::Duct)?;
        match &self.branch_coupling {
            Some(coupling) => coupling.response(main, &self.branches),
            None => Ok(main),
        }
    }

    // The caller owns all private candidates; no step below publishes state.
    pub(super) fn branch_inputs(
        &self,
        external: &[PortInput],
    ) -> Result<(Vec<PortInput>, Vec<f64>), AcousticRealizeError> {
        let Some(coupling) = &self.branch_coupling else {
            return Ok((external.to_vec(), Vec::new()));
        };
        let main = self.propagator.multiport_response(&self.state)
            .map_err(AcousticRealizeError::Duct)?;
        let branches = coupling::responses(&self.branches)?;
        let flows = coupling.flows(&main, &branches, external)?;
        let inputs = main
            .ports
            .iter()
            .map(|&id| {
                if let Some(index) = self.branches.iter().position(|branch| branch.port == id) {
                    Ok(PortInput {
                        id,
                        effort: -flows[index],
                    })
                } else {
                    external
                        .iter()
                        .find(|input| input.id == id)
                        .copied()
                        .ok_or_else(|| invalid("DG bore has no owned exterior channel"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((inputs, flows))
    }

    pub(super) fn step_branches(
        &self,
        external: &[PortInput],
        main_inputs: &[PortInput],
        junction_flows: &[f64],
        output: &mut DgStepOutput,
    ) -> Result<(Vec<DgState>, Vec<DgBranchFrame>, Vec<f64>), AcousticRealizeError> {
        if self.branches.is_empty() { return Ok((Vec::new(), Vec::new(), Vec::new())); }
        let main_response = self.propagator.multiport_response(&self.state)
            .map_err(AcousticRealizeError::Duct)?;
        let mut original_pressure_scales = main_response.ports.iter().map(|&port|
            pressure_scale(&main_response, port, main_inputs)).collect::<Result<Vec<_>, _>>()?;
        let mut states = Vec::with_capacity(self.branches.len());
        let mut frames = Vec::with_capacity(self.branches.len());
        for (branch, &flow) in self.branches.iter().zip(junction_flows) {
            let mut inputs = vec![PortInput { id: branch.inlet, effort: flow }];
            if let Some(outlet) = branch.outlet {
                let distal = external.iter().find(|input| input.id == branch.port)
                    .ok_or_else(|| invalid("DG branch exterior flow was omitted"))?;
                inputs.push(PortInput { id: outlet, effort: distal.effort });
            }
            let response = branch.propagator.multiport_response(&branch.state)
                .map_err(AcousticRealizeError::Duct)?;
            let mut state = branch.state.clone();
            let step = branch.propagator.step(&mut state, &inputs)
                .map_err(AcousticRealizeError::Duct)?;
            branch.propagator.sync_nodal(&mut state).map_err(AcousticRealizeError::Duct)?;
            validate_candidate(&state, &step, inputs.len())?;
            if step.ledger.step != output.ledger.step {
                return Err(invalid("DG branch and bore clocks diverged"));
            }
            let bore_port = *output.port_outputs.iter().find(|port| port.id == branch.port)
                .ok_or_else(|| invalid("DG bore omitted a branch junction"))?;
            let branch_inlet = *step.port_outputs.iter().find(|port| port.id == branch.inlet)
                .ok_or_else(|| invalid("DG branch omitted its inlet work"))?;
            let junction_pressure_scale = pressure_scale(&main_response, branch.port, main_inputs)?
                + pressure_scale(&response, branch.inlet, &inputs)?;
            let bore_work_j = flow * bore_port.delta_x;
            let branch_work_j = flow * branch_inlet.delta_x;
            let residual = bore_work_j - branch_work_j;
            let work_scale = junction_pressure_scale * flow.abs() * self.propagator.dt();
            if ![junction_pressure_scale, bore_work_j, branch_work_j, residual, work_scale]
                .iter().all(|value| value.is_finite())
                || (bore_port.conjugate_rate - branch_inlet.conjugate_rate).abs()
                    > 1e-10 * junction_pressure_scale.max(1e-300)
                || residual.abs() > 1e-10 * work_scale.max(1e-300)
            {
                return Err(invalid("DG branch junction disagrees on same-step pressure or work"));
            }
            let external_slot = main_response.ports.iter().position(|id| *id == branch.port)
                .ok_or_else(|| invalid("DG branch lost its original external pressure basis"))?;
            if let Some(outlet) = branch.outlet {
                original_pressure_scales[external_slot] = pressure_scale(&response, outlet, &inputs)?;
                let distal_port = *step.port_outputs.iter().find(|port| port.id == outlet)
                    .ok_or_else(|| invalid("DG branch omitted its exterior output"))?;
                let port = output.port_outputs.iter_mut().find(|port| port.id == branch.port)
                    .ok_or_else(|| invalid("DG branch lost its external channel identity"))?;
                *port = PortOutput { id: branch.port, ..distal_port };
            } else {
                // The inlet exchange remains in DgBranchFrame; a rigid distal
                // wall contributes no fictitious exterior pressure/flow row.
                output.port_outputs.retain(|port| port.id != branch.port);
            }
            let ledger = &step.ledger;
            output.ledger.stored_energy_j += ledger.stored_energy_j;
            output.ledger.port_work_step_j += ledger.port_work_step_j;
            output.ledger.cumulative_port_work_j += ledger.cumulative_port_work_j;
            output.ledger.cumulative_dissipated_j += ledger.cumulative_dissipated_j;
            output.ledger.cumulative_parameter_work_j += ledger.cumulative_parameter_work_j;
            output.ledger.energy_balance_residual_j += ledger.energy_balance_residual_j;
            output.ledger.balance_checked &= ledger.balance_checked;
            frames.push(DgBranchFrame {
                port: branch.port, junction_flow_m3_s: flow, bore_port, branch_inlet,
                bore_work_j, branch_work_j, junction_work_residual_j: residual, duct: step,
            });
            states.push(state);
        }
        // Every participant was constructed cold. Retain the signed original
        // residual sum; internal cancellation roundoff is not clipped away.
        output.ledger.relative_energy_drift = output.ledger.energy_balance_residual_j.abs()
            / (output.ledger.cumulative_port_work_j.abs()
                + output.ledger.cumulative_parameter_work_j.abs()).max(1e-300);
        let original_pressure_scales = main_response
            .ports
            .iter()
            .zip(original_pressure_scales)
            .filter_map(|(port, scale)| {
                (!self
                    .branches
                    .iter()
                    .any(|branch| branch.port == *port && branch.outlet.is_none()))
                .then_some(scale)
            })
            .collect();
        Ok((states, frames, original_pressure_scales))
    }
}

pub(super) fn admit_sources(
    main: &DgDuct,
    branches: &[DgBranch],
    dt: f64,
    maximum_dofs: usize,
) -> Result<(), AcousticRealizeError> {
    let mut total = main.total_dofs();
    for (index, branch) in branches.iter().enumerate() {
        if branches[..index].iter().any(|previous| previous.port == branch.port)
            || !main.side_ports().iter().any(|side| side.id == branch.port)
            || !branch.duct.side_ports().is_empty()
            || branch.duct.inlet() == branch.duct.outlet()
            || !matches!((branch.duct.inlet(), branch.duct.outlet()),
                (DgBoundary::Port(_), DgBoundary::Port(_) | DgBoundary::Closed))
        {
            return Err(invalid("DG branches require distinct original side junctions and a port or closed distal boundary"));
        }
        let source = main.side_port_source(branch.port)
            .map_err(|_| invalid("DG branch junction has no original local material source"))?;
        if source.gas.pressure.to_bits() != branch.duct.inlet_gas().pressure.to_bits() {
            return Err(invalid("DG branch junction requires equal equilibrium material pressure"));
        }
        if branch.duct.viscothermal_band_rad_s().is_some_and(|band|
            band[1] * dt > core::f64::consts::PI)
        {
            return Err(invalid("DG branch viscothermal band exceeds the mechanical Nyquist frequency"));
        }
        total = total.checked_add(branch.duct.total_dofs())
            .ok_or_else(|| invalid("DG branch state allowance overflow"))?;
    }
    if total > maximum_dofs {
        return Err(invalid("DG bore plus branches acoustic and memory states exceed the DOF budget"));
    }
    Ok(())
}

fn pressure_scale(
    response: &DgMultiportResponse, port: PortId, inputs: &[PortInput],
) -> Result<f64, AcousticRealizeError> {
    let n = response.ports.len();
    let row = response.ports.iter().position(|id| *id == port)
        .ok_or_else(|| invalid("DG pressure audit lost its original port"))?;
    let mut scale = response.history_pressure_pa[row].abs();
    for input in inputs {
        let col = response.ports.iter().position(|id| *id == input.id)
            .ok_or_else(|| invalid("DG pressure audit received a foreign port"))?;
        scale += (response.impedance_pa_s_m3[row * n + col] * input.effort).abs();
    }
    Ok(scale)
}

pub(super) fn validate_candidate(
    state: &DgState, output: &DgStepOutput, ports: usize,
) -> Result<(), AcousticRealizeError> {
    let ledger = &output.ledger;
    if output.port_outputs.len() != ports
        || output.port_outputs.iter().any(|port|
            !port.conjugate_rate.is_finite() || !port.delta_x.is_finite()
                || !port.end_conjugate.is_finite())
        || ![
            state.stored_energy_j(), ledger.stored_energy_j, ledger.port_work_step_j,
            ledger.cumulative_port_work_j, ledger.cumulative_dissipated_j,
            ledger.energy_balance_residual_j, ledger.relative_energy_drift,
        ].iter().all(|value| value.is_finite())
        || state.pressure().iter().chain(state.volume_flow()).any(|value| !value.is_finite())
    {
        return Err(invalid("DG port candidate left the finite pressure/state/work domain"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
