//! Source-bound exact mean-port loads shared by nonlinear acoustic couplers.
//!
//! The original DG owner retains geometry, material gas, wall memory, storage
//! and propagation. An optional admitted rational port group exchanges same-step
//! mean pressures and opposite volume flows through its original MIMO ZOH owner.
//! Candidates remain private until every participant accepts the shared work.

/// Authored geometry, material gas and numerical declarations for this load.
pub mod input;

/// Computed exterior radiation from the actual endpoint geometry and gas.
pub mod radiation;

mod terminal;
mod continuous;
pub(crate) use continuous::DgTerminalTrajectory;
use terminal::{TerminalCoupling, TerminalRelation};
mod mouth;
pub use mouth::DgDrivenMouthFrame;

mod branch;
pub use branch::{DgBranch, DgBranchFrame, DgBranchLoad};
use branch::BranchCoupling;

mod channel;
pub use channel::{
    DgTerminalChannelFrame, DgTerminalChannelInput, DgTerminalChannelKind,
    DgTerminalChannelLoad, MAX_DG_TERMINAL_TRANSVERSE_CELLS,
};

use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledRecord;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgMultiportResponse, DgPortResponse, DgPropagator, DgState,
    DgStepOutput, PortId, PortInput,
};
use fs_vfit::mimo::load::PreparedMimoLoad;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Accepted exchange from an original DG port into its rational load.
/// Work is SIGNED exchange, not a dissipation or rational-storage certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgTerminalFrame {
    /// Stable external channel ID. For a finite side branch this denotes its
    /// actual outer outlet; its original bore-side ID is the internal junction.
    pub port: PortId,
    /// Interval-mean outward flow [m³/s]; held only on the mean-port backend.
    pub outward_flow_m3_s: f64,
    /// Mean outlet pressure from the original DG step [Pa].
    pub mean_pressure_pa: f64,
    /// Actual integral of DG pressure times outward flow [J].
    pub port_work_j: f64,
    /// Actual integral of original rational pressure times its input flow [J].
    pub load_work_j: f64,
    /// Work into the optional original finite gap channel [J]. This contains
    /// kinetic storage and is not all heat. Zero on a direct or closed port.
    pub channel_work_j: f64,
    /// DG work minus original channel work and original exterior work [J].
    pub port_work_residual_j: f64,
}

// Preserve a scalar exchange's original value and sum groups in declared order.
pub(crate) fn sum_terminal_work(frames: &[DgTerminalFrame]) -> f64 {
    frames.split_first().map_or(0.0, |(first, rest)|
        rest.iter().fold(first.port_work_j - first.channel_work_j,
            |sum, frame| sum + (frame.port_work_j - frame.channel_work_j)))
}

/// One declared inlet, one original propagator and one accepted physical history.
/// Read its effective mean impedance before constructing a nonlinear source.
#[derive(Debug, Clone)]
pub struct DgPortLoad {
    duct: DgDuct,
    propagator: DgPropagator,
    state: DgState,
    inlet: PortId,
    impedance_pa_s_m3: f64,
    terminal: Option<PreparedMimoLoad>,
    terminal_coupling: Option<TerminalCoupling>,
    last_terminal_frame: Option<DgTerminalFrame>,
    last_terminal_frames: Vec<DgTerminalFrame>,
    last_driven_mouth_frame: Option<DgDrivenMouthFrame>,
    radiation_record: Option<ComputedBaffledRecord>,
    branches: Vec<DgBranchLoad>,
    branch_coupling: Option<BranchCoupling>,
    last_branch_frames: Vec<DgBranchFrame>,
    terminal_channels: Vec<DgTerminalChannelLoad>,
    last_terminal_channel_frames: Vec<DgTerminalChannelFrame>,
}

impl DgPortLoad {
    /// Prepare an initially quiescent duct with one inlet port and a closed,
    /// ideal-open, or matched outlet. The DOF allowance covers acoustic AND wall
    /// states before construction of the original dense exponential.
    ///
    /// The complete material wall band must lie below mechanical Nyquist.
    /// A zero return-period mean impedance is refused, never clamped.
    pub fn new(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
    ) -> Result<Self, AcousticRealizeError> {
        Self::prepare(duct, time_step_s, exponential_tolerance, maximum_dofs, Vec::new(), None, Vec::new())
    }

    /// Attach one admitted causal impedance to the declared outlet port.
    /// The source receives the Schur complement of the full DG mean-port matrix
    /// and the original rational mean law; no feedback sample is delayed.
    ///
    /// The terminal's own admission checks shape, finite/stable/proper source,
    /// full mechanical Nyquist bandwidth, state/memory budgets, continuous
    /// passivity and actual sampled-mean passivity. Its original report classes
    /// remain available through terminal(); sampled admission is still GridOnly.
    ///
    /// Each participant integrates exactly for its held boundary flow. This
    /// work-conjugate mean coupling is not the joint continuous exponential;
    /// interconnection-clock refinement remains a separate accuracy requirement.
    ///
    /// # Errors
    /// Missing/distinct outlet, multiple terminal channels, clock mismatch,
    /// nonquiescent terminal memory, or nonpositive effective mean impedance,
    /// in addition to the original DG geometry/material/exponential admissions.
    pub fn new_with_terminal(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        terminal: PreparedMimoLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(outlet) = duct.outlet() else {
            return Err(invalid("DG terminal requires a declared axial outlet port"));
        };
        Self::new_with_terminal_ports(
            duct, time_step_s, exponential_tolerance, maximum_dofs, vec![outlet], terminal,
        )
    }

    /// Attach one original MIMO impedance to an explicit ordered terminal group.
    /// Every declared DG port except the source inlet must appear exactly once;
    /// no omitted side port or outlet silently receives zero flow. Load rows,
    /// columns, pole histories and returned frames use this exact channel order.
    ///
    /// The complete same-step relation is `(Zdd + L) u = Hd - Hl + Zd0 Qin`,
    /// where `u` flows OUT of DG and into the load. Cross terms are preserved.
    /// The scalar case keeps the original divide arithmetic. Larger constant
    /// systems use the existing prepared LU once at construction, with a finite
    /// infinity-norm condition bound of 1e10; runtime histories use that solve.
    ///
    /// Each original owner integrates exactly under held flows. As with the
    /// scalar outlet, this is a work-conjugate mean connection, not a claim of
    /// the joint continuous exponential. Refine the coupling clock separately.
    pub fn new_with_terminal_ports(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        ports: Vec<PortId>,
        terminal: PreparedMimoLoad,
    ) -> Result<Self, AcousticRealizeError> {
        Self::prepare(
            duct, time_step_s, exponential_tolerance, maximum_dofs, ports, Some(terminal), Vec::new(),
        )
    }

    fn prepare(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        terminal_ports: Vec<PortId>,
        terminal: Option<PreparedMimoLoad>,
        branches: Vec<DgBranch>,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(inlet) = duct.inlet() else {
            return Err(invalid("DG load requires one declared inlet port"));
        };
        if let Some(load) = &terminal {
            if load.time_step_s().to_bits() != time_step_s.to_bits()
                || load.pole_states().iter().any(|z| z.re != 0.0 || z.im != 0.0)
            {
                return Err(invalid("DG terminal requires quiescent original memory on the same clock"));
            }
        }
        if !branches.is_empty() {
            branch::admit_sources(&duct, &branches, time_step_s, maximum_dofs)?;
        }
        if duct.total_dofs() > maximum_dofs {
            return Err(invalid("DG acoustic plus memory states exceed the declared DOF budget"));
        }
        if let Some(band) = duct.viscothermal_band_rad_s() {
            if band[1] * time_step_s > core::f64::consts::PI {
                return Err(invalid("DG viscothermal band exceeds the mechanical Nyquist frequency"));
            }
        }
        let propagator = duct
            .propagator(time_step_s, exponential_tolerance)
            .map_err(AcousticRealizeError::Duct)?;
        let state = duct.zero_state();
        let response = propagator.multiport_response(&state)
            .map_err(AcousticRealizeError::Duct)?;
        let branches: Vec<_> = branches.into_iter().map(|source|
            DgBranchLoad::prepare(source, time_step_s, exponential_tolerance)).collect::<Result<_, _>>()?;
        let branch_coupling = if branches.is_empty() { None }
            else { Some(BranchCoupling::new(&response, &branches)?) };
        let response = match &branch_coupling {
            Some(coupling) => coupling.response(response, &branches)?,
            None => response,
        };
        let terminal_coupling = match &terminal {
            Some(load) => Some(TerminalCoupling::new(&response, inlet, terminal_ports, load)?),
            None => {
                if response.ports != [inlet] {
                    return Err(invalid("DG load requires an explicit owner for every non-inlet flow port"));
                }
                None
            }
        };
        let mut load = Self {
            duct, propagator, state, inlet, impedance_pa_s_m3: 0.0,
            terminal, terminal_coupling, last_terminal_frame: None,
            last_terminal_frames: Vec::new(), last_driven_mouth_frame: None, radiation_record: None,
            branches, branch_coupling, last_branch_frames: Vec::new(),
            terminal_channels: Vec::new(), last_terminal_channel_frames: Vec::new(),
        };
        let impedance = load.port_response()?.impedance_pa_s_m3;
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid("DG load requires strictly positive exact mean-port impedance"));
        }
        load.impedance_pa_s_m3 = impedance;
        Ok(load)
    }

    /// Effective mean source effort per held inlet flow, including all reaction.
    /// A driven mouth uses the physical inlet-minus-mouth pressure difference.
    #[must_use]
    pub const fn impedance_pa_s_m3(&self) -> f64 { self.impedance_pa_s_m3 }

    /// Original geometry and material source.
    #[must_use]
    pub const fn duct(&self) -> &DgDuct { &self.duct }

    /// Original exponential, admitted clock and certificate.
    #[must_use]
    pub const fn propagator(&self) -> &DgPropagator { &self.propagator }

    /// Accepted MAIN bore state, including current physical nodal fields.
    /// Finite branches retain independent original states through branches().
    #[must_use]
    pub const fn state(&self) -> &DgState { &self.state }

    /// The unique inlet port receiving positive flow into the duct.
    #[must_use]
    pub const fn inlet(&self) -> PortId { self.inlet }

    /// Original admitted terminal, source matrices, pole history and reports.
    #[must_use]
    pub fn terminal(&self) -> Option<&PreparedMimoLoad> { self.terminal.as_ref() }

    /// Bound DG terminal ports in their original relative MIMO channel order.
    /// An optional driven exterior mouth occupies the separate leading channel.
    #[must_use]
    pub fn terminal_ports(&self) -> &[PortId] {
        self.terminal_coupling.as_ref().map_or(&[], |coupling| coupling.ports.as_slice())
    }

    /// Physical geometry, exterior material and original computed-load reports.
    /// Generic authored terminals do not manufacture a radiation source record.
    #[must_use]
    pub fn radiation_record(&self) -> Option<&ComputedBaffledRecord> {
        self.radiation_record.as_ref()
    }

    /// Last accepted AXIAL outlet exchange; absent before the first step or
    /// when that actual outlet is not loaded. Side-only groups do not pretend
    /// that their first channel is the axial outlet.
    #[must_use]
    pub const fn last_terminal_frame(&self) -> Option<&DgTerminalFrame> {
        self.last_terminal_frame.as_ref()
    }

    /// All accepted exchanges in terminal_ports order. Empty before the first
    /// accepted step; a failed or abandoned trial preserves the complete slice.
    #[must_use]
    pub fn last_terminal_frames(&self) -> &[DgTerminalFrame] { &self.last_terminal_frames }

    fn terminal_relation(
        &self, terminal: &PreparedMimoLoad,
    ) -> Result<TerminalRelation, AcousticRealizeError> {
        let response = self.multiport_response()?;
        self.terminal_coupling.as_ref()
            .ok_or_else(|| invalid("DG terminal lost its admitted port group"))?
            .relation(response, terminal, &self.terminal_channels)
    }

    /// Affine source effort from all accepted acoustic, wall and load histories.
    /// This is inlet pressure for ordinary terminals. An explicitly attached
    /// driven exterior mouth instead returns inlet minus exterior-mouth pressure.
    /// No participant advances while this relation is inspected.
    pub fn port_response(&self) -> Result<DgPortResponse, AcousticRealizeError> {
        match &self.terminal {
            Some(terminal) => Ok(self.terminal_relation(terminal)?.inlet),
            None => self.propagator.port_response(&self.state, self.inlet)
                .map_err(AcousticRealizeError::Duct),
        }
    }

    /// Retain private DG and terminal candidates until the nonlinear participant
    /// validates its own state and work. Dropping this trial publishes nothing.
    pub(crate) fn preview_step(
        &mut self, volume_flow_m3_s: f64,
    ) -> Result<DgPortTrial<'_>, AcousticRealizeError> {
        if !volume_flow_m3_s.is_finite() {
            return Err(invalid("DG inlet flow must be finite"));
        }
        let relation = self.terminal.as_ref()
            .map(|terminal| self.terminal_relation(terminal)).transpose()?;
        let mut terminal = self.terminal.clone();
        let mut inputs = vec![PortInput { id: self.inlet, effort: volume_flow_m3_s }];
        let count = self.terminal_ports().len();
        let mut outward_flows = vec![0.0; count];
        let load_offset = usize::from(self.terminal_coupling.as_ref()
            .is_some_and(|coupling| coupling.driven_mouth));
        let mut load_pressure = vec![0.0; count + load_offset];
        if let (Some(relation), Some(candidate), Some(coupling)) =
            (&relation, &mut terminal, &self.terminal_coupling)
        {
            for (i, outward) in outward_flows.iter_mut().enumerate() {
                *outward = relation.outward_history_m3_s[i]
                    + coupling.outward_per_inlet[i] * volume_flow_m3_s;
            }
            let mut load_flows = Vec::with_capacity(count + load_offset);
            if load_offset != 0 { load_flows.push(-volume_flow_m3_s); }
            load_flows.extend_from_slice(&outward_flows);
            candidate.step_mean_into(&load_flows, &mut load_pressure)
                .map_err(|_| invalid("DG rational terminal candidate left the finite pressure/state domain"))?;
            for (id, outward) in coupling.ports.iter().zip(&outward_flows) {
                inputs.push(PortInput { id: *id, effort: -*outward });
            }
        }
        let (main_inputs, junction_flows) = self.branch_inputs(&inputs)?;
        let mut state = self.state.clone();
        let mut output = self.propagator.step(&mut state, &main_inputs)
            .map_err(AcousticRealizeError::Duct)?;
        self.propagator.sync_nodal(&mut state).map_err(AcousticRealizeError::Duct)?;
        branch::validate_candidate(&state, &output, inputs.len())?;
        let (branch_states, branch_frames, original_pressure_scales) =
            self.step_branches(&inputs, &main_inputs, &junction_flows, &mut output)?;
        if !branch_frames.is_empty() { branch::validate_candidate(&state, &output, inputs.len())?; }
        let inlet_port = output.port_outputs.iter().find(|port| port.id == self.inlet)
            .ok_or_else(|| invalid("DG candidate omitted its original inlet"))?;
        let inlet_work_j = volume_flow_m3_s * inlet_port.delta_x;
        let driven_mouth_frame = if load_offset != 0 {
            let source = &relation.as_ref().expect("driven mouth requires its admitted terminal relation").inlet;
            Some(DgDrivenMouthFrame::new(volume_flow_m3_s, load_pressure[0],
                inlet_port.conjugate_rate, inlet_work_j, source, self.propagator.dt())?)
        } else { None };
        let source_work_j = driven_mouth_frame.as_ref()
            .map_or(inlet_work_j, |frame| frame.source_work_j);
        let (terminal_channels, terminal_channel_frames) = self.step_terminal_channels(
            &output, &outward_flows, &load_pressure, load_offset)?;
        let mut terminal_frames = Vec::with_capacity(count);
        if let (Some(relation), Some(coupling)) = (relation, &self.terminal_coupling) {
            let h = &relation.response.history_pressure_pa;
            let z = &relation.response.impedance_pa_s_m3;
            let width = relation.response.ports.len();
            for (i, &id) in coupling.ports.iter().enumerate() {
                let port = output.port_outputs.iter().find(|port| port.id == id)
                    .ok_or_else(|| invalid("DG candidate omitted an owned terminal port"))?;
                let row = coupling.indices[i];
                let mut pressure_scale = h[row].abs()
                    + (z[width * row + coupling.inlet_index] * volume_flow_m3_s).abs();
                for (j, &column) in coupling.indices.iter().enumerate() {
                    pressure_scale += (z[width * row + column] * outward_flows[j]).abs();
                }
                let terminal_pressure = load_pressure[i + load_offset];
                pressure_scale += terminal_pressure.abs();
                let channel = terminal_channel_frames.iter().find(|frame| frame.port == id);
                let channel_pressure = channel.map_or(0.0, |frame| frame.constitutive_pressure_drop_pa);
                let channel_work_j = channel.map_or(0.0, |frame| frame.pressure_work_j);
                let closed = channel.is_some_and(|frame| frame.step.is_none());
                pressure_scale += channel_pressure.abs();
                // Branch condensation can cancel large internal terms. Keep the
                // original owner terms in the same fixed relative work audit.
                if let Some(&original_scale) = original_pressure_scales.get(row) {
                    pressure_scale += original_scale;
                }
                let outward_flow = outward_flows[i];
                let port_work_j = outward_flow * port.delta_x;
                let load_work_j = terminal_pressure * (outward_flow * self.propagator.dt());
                let port_work_residual_j = (port_work_j - channel_work_j) - load_work_j;
                let work_scale = pressure_scale * outward_flow.abs() * self.propagator.dt()
                    + channel.map_or(0.0, DgTerminalChannelFrame::work_scale_j);
                if ![
                    outward_flow, terminal_pressure, pressure_scale, port_work_j,
                    load_work_j, channel_work_j, port_work_residual_j, work_scale,
                ].iter().all(|value| value.is_finite())
                    || (!closed && ((port.conjugate_rate - terminal_pressure) - channel_pressure).abs()
                        > 1e-10 * pressure_scale.max(1e-300))
                    || port_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
                {
                    return Err(invalid("DG and rational terminal disagree on same-step pressure or work"));
                }
                terminal_frames.push(DgTerminalFrame {
                    port: id, outward_flow_m3_s: outward_flow,
                    mean_pressure_pa: port.conjugate_rate,
                    port_work_j, load_work_j, channel_work_j, port_work_residual_j,
                });
            }
        }
        let terminal_frame = match self.duct.outlet() {
            DgBoundary::Port(outlet) => terminal_frames.iter().find(|frame| frame.port == outlet).copied(),
            _ => None,
        };
        let terminal_work = sum_terminal_work(&terminal_frames);
        channel::include_channel_ledgers(&mut output, &terminal_channels)?;
        let ledger = &output.ledger;
        let balance = ledger.port_work_step_j - (inlet_work_j - terminal_work);
        let terminal_work_scale = terminal_frames.iter().map(|frame| frame.port_work_j.abs()).sum::<f64>();
        let mut scale = ledger.port_work_step_j.abs() + inlet_work_j.abs() + terminal_work_scale;
        for frame in &branch_frames {
            scale += frame.bore_work_j.abs() + frame.branch_work_j.abs();
        }
        for frame in &terminal_channel_frames { scale += frame.work_scale_j(); }
        if ![inlet_work_j, balance, scale].iter().all(|value| value.is_finite())
            || balance.abs() > 1e-10 * scale.max(1e-300)
        {
            return Err(invalid("aggregate DG/channel port work does not equal inlet minus exterior work"));
        }
        Ok(DgPortTrial { owner: self, state, output, terminal, terminal_frame, terminal_frames,
            inlet_work_j, source_work_j, driven_mouth_frame, branch_states, branch_frames,
            terminal_channels, terminal_channel_frames })
    }
}

/// A trial borrows the only accepted histories. Dropping it publishes nothing.
pub(crate) struct DgPortTrial<'a> {
    owner: &'a mut DgPortLoad,
    state: DgState,
    output: DgStepOutput,
    terminal: Option<PreparedMimoLoad>,
    terminal_frame: Option<DgTerminalFrame>,
    terminal_frames: Vec<DgTerminalFrame>,
    inlet_work_j: f64,
    source_work_j: f64,
    driven_mouth_frame: Option<DgDrivenMouthFrame>,
    branch_states: Vec<DgState>,
    branch_frames: Vec<DgBranchFrame>,
    terminal_channels: Vec<DgTerminalChannelLoad>,
    terminal_channel_frames: Vec<DgTerminalChannelFrame>,
}

impl DgPortTrial<'_> {
    pub(crate) const fn previous_state(&self) -> &DgState { &self.owner.state }
    pub(crate) const fn state(&self) -> &DgState { &self.state }
    /// Actual external DG pressure observations and the aggregate bore, branch
    /// and gap ledger. With finite terminal gaps, these are the duct-side
    /// pressures; terminal/channel frames retain the exterior pressure/work.
    /// Internal junction work remains visible in branch_frames().
    pub(crate) const fn output(&self) -> &DgStepOutput { &self.output }
    pub(crate) fn previous_stored_energy_j(&self) -> f64 { self.owner.stored_energy_j() }
    pub(crate) fn previous_cumulative_dissipated_j(&self) -> f64 {
        self.owner.cumulative_dissipated_j()
    }
    pub(crate) fn branch_frames(&self) -> &[DgBranchFrame] { &self.branch_frames }
    pub(crate) const fn inlet_work_j(&self) -> f64 { self.inlet_work_j }
    pub(crate) const fn source_work_j(&self) -> f64 { self.source_work_j }
    pub(crate) const fn driven_mouth_frame(&self) -> Option<&DgDrivenMouthFrame> {
        self.driven_mouth_frame.as_ref()
    }
    pub(crate) const fn terminal_frame(&self) -> Option<&DgTerminalFrame> {
        self.terminal_frame.as_ref()
    }
    pub(crate) fn terminal_frames(&self) -> &[DgTerminalFrame] { &self.terminal_frames }
    pub(crate) fn terminal_channel_frames(&self) -> &[DgTerminalChannelFrame] {
        &self.terminal_channel_frames
    }

    /// Publish every already validated candidate; all fallible checks precede this.
    pub(crate) fn accept(self) -> DgStepOutput {
        self.owner.state = self.state;
        self.owner.terminal = self.terminal;
        self.owner.last_terminal_frame = self.terminal_frame;
        self.owner.last_terminal_frames = self.terminal_frames;
        self.owner.last_driven_mouth_frame = self.driven_mouth_frame;
        for (branch, state) in self.owner.branches.iter_mut().zip(self.branch_states) {
            branch.state = state;
        }
        self.owner.last_branch_frames = self.branch_frames;
        self.owner.terminal_channels = self.terminal_channels;
        self.owner.last_terminal_channel_frames = self.terminal_channel_frames;
        self.output
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod terminal_tests;
#[cfg(test)]
mod terminal_group_tests;

