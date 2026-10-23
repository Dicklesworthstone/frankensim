//! Source-bound exact mean-port loads shared by nonlinear acoustic couplers.
//!
//! The original DG owner retains geometry, material gas, wall memory, storage
//! and propagation. An optional admitted rational outlet exchanges same-step
//! mean pressure and opposite volume flow through its original MIMO ZOH owner.
//! Candidates remain private until every participant accepts the shared work.

/// Authored geometry, material gas and numerical declarations for this load.
pub mod input;

use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgMultiportResponse, DgPortResponse, DgPropagator, DgState,
    DgStepOutput, PortId, PortInput,
};
use fs_vfit::mimo::load::PreparedMimoLoad;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Same-step exchange from the DG outlet into its original rational load.
/// Work is SIGNED exchange, not a dissipation or rational-storage certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgTerminalFrame {
    /// The original outlet port; its DG input has the opposite flow sign.
    pub port: PortId,
    /// Held volume flow out of the duct and into the terminal [m³/s].
    pub outward_flow_m3_s: f64,
    /// Mean outlet pressure from the original DG step [Pa].
    pub mean_pressure_pa: f64,
    /// Exact DG outlet impulse times outward flow [J].
    pub port_work_j: f64,
    /// Original rational mean pressure times held flow and duration [J].
    pub load_work_j: f64,
    /// Difference of the independently evaluated port works [J].
    pub port_work_residual_j: f64,
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
    last_terminal_frame: Option<DgTerminalFrame>,
}

// This is only the algebraic elimination of the outlet's same-step mean law.
// The two existing owners still advance their complete original coordinates.
struct TerminalRelation {
    inlet: DgPortResponse,
    outlet: PortId,
    outward_history_m3_s: f64,
    outward_per_inlet: f64,
    response: DgMultiportResponse,
    inlet_index: usize,
    outlet_index: usize,
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
        Self::prepare(duct, time_step_s, exponential_tolerance, maximum_dofs, None)
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
        Self::prepare(
            duct, time_step_s, exponential_tolerance, maximum_dofs, Some(terminal),
        )
    }

    fn prepare(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        terminal: Option<PreparedMimoLoad>,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(inlet) = duct.inlet() else {
            return Err(invalid("DG load requires one declared inlet port"));
        };
        match (&terminal, duct.outlet()) {
            (None, DgBoundary::Port(_)) => {
                return Err(invalid("DG load requires an owner for its second unknown port"));
            }
            (Some(load), DgBoundary::Port(outlet))
                if outlet != inlet && load.ports() == 1
                    && load.time_step_s().to_bits() == time_step_s.to_bits()
                    && load.pole_states().iter().all(|z| z.re == 0.0 && z.im == 0.0) => {}
            (Some(_), _) => {
                return Err(invalid("DG terminal requires a distinct outlet and one quiescent channel on the same clock"));
            }
            (None, _) => {}
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
        let mut load = Self {
            duct, propagator, state, inlet, impedance_pa_s_m3: 0.0,
            terminal, last_terminal_frame: None,
        };
        let impedance = load.port_response()?.impedance_pa_s_m3;
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid("DG load requires strictly positive exact mean-port impedance"));
        }
        load.impedance_pa_s_m3 = impedance;
        Ok(load)
    }

    /// Effective mean pressure per held inlet flow, including outlet reaction.
    #[must_use]
    pub const fn impedance_pa_s_m3(&self) -> f64 { self.impedance_pa_s_m3 }

    /// Original geometry and material source.
    #[must_use]
    pub const fn duct(&self) -> &DgDuct { &self.duct }

    /// Original exponential, admitted clock and certificate.
    #[must_use]
    pub const fn propagator(&self) -> &DgPropagator { &self.propagator }

    /// Complete accepted DG state, including current physical nodal fields.
    #[must_use]
    pub const fn state(&self) -> &DgState { &self.state }

    /// The unique inlet port receiving positive flow into the duct.
    #[must_use]
    pub const fn inlet(&self) -> PortId { self.inlet }

    /// Original admitted terminal, source matrices, pole history and reports.
    #[must_use]
    pub fn terminal(&self) -> Option<&PreparedMimoLoad> { self.terminal.as_ref() }

    /// Last accepted outlet exchange; absent before the first step or without
    /// an explicit rational outlet. This flow is suitable for its bound observer.
    #[must_use]
    pub const fn last_terminal_frame(&self) -> Option<&DgTerminalFrame> {
        self.last_terminal_frame.as_ref()
    }

    fn terminal_relation(
        &self, terminal: &PreparedMimoLoad,
    ) -> Result<TerminalRelation, AcousticRealizeError> {
        let DgBoundary::Port(outlet) = self.duct.outlet() else {
            return Err(invalid("DG terminal lost its declared outlet"));
        };
        let response = self.propagator.multiport_response(&self.state)
            .map_err(AcousticRealizeError::Duct)?;
        if response.ports.len() != 2 {
            return Err(invalid("DG terminal needs exactly its original inlet and outlet"));
        }
        let inlet_index = response.ports.iter().position(|id| *id == self.inlet)
            .ok_or_else(|| invalid("DG inlet missing from its original port basis"))?;
        let outlet_index = response.ports.iter().position(|id| *id == outlet)
            .ok_or_else(|| invalid("DG outlet missing from its original port basis"))?;
        let z = &response.impedance_pa_s_m3;
        let h = &response.history_pressure_pa;
        let mut load_history = [0.0];
        terminal.mean_history_output_into(&mut load_history)
            .map_err(|_| invalid("DG terminal mean history is outside its admitted finite domain"))?;
        let denominator = z[2 * outlet_index + outlet_index] + terminal.mean_direct_term()[0];
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(invalid("DG terminal mean solve requires positive finite combined outlet impedance"));
        }
        // Both DG inputs are inward; the terminal input is outward.
        // u = (Ho - Hl + Zoi Qi)/(Zoo + Zl)
        // Pi = Hi + Zii Qi - Zio u. Cross terms need not be symmetric.
        let outward_history_m3_s = (h[outlet_index] - load_history[0]) / denominator;
        let outward_per_inlet = z[2 * outlet_index + inlet_index] / denominator;
        let inlet = DgPortResponse {
            history_pressure_pa: h[inlet_index]
                - z[2 * inlet_index + outlet_index] * outward_history_m3_s,
            impedance_pa_s_m3: z[2 * inlet_index + inlet_index]
                - z[2 * inlet_index + outlet_index] * outward_per_inlet,
        };
        if ![
            outward_history_m3_s, outward_per_inlet,
            inlet.history_pressure_pa, inlet.impedance_pa_s_m3,
        ].iter().all(|value| value.is_finite()) || inlet.impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("DG terminal elimination left the finite positive inlet-load domain"));
        }
        Ok(TerminalRelation {
            inlet, outlet, outward_history_m3_s, outward_per_inlet,
            response, inlet_index, outlet_index,
        })
    }

    /// Affine inlet pressure from all accepted acoustic, wall and terminal
    /// histories. No participant advances while this relation is inspected.
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
        let mut outward_flow = 0.0;
        let mut load_pressure = [0.0];
        if let (Some(relation), Some(candidate)) = (&relation, &mut terminal) {
            outward_flow = relation.outward_history_m3_s
                + relation.outward_per_inlet * volume_flow_m3_s;
            candidate.step_mean_into(&[outward_flow], &mut load_pressure)
                .map_err(|_| invalid("DG rational outlet candidate left the finite pressure/state domain"))?;
            inputs.push(PortInput { id: relation.outlet, effort: -outward_flow });
        }
        let mut state = self.state.clone();
        let output = self.propagator.step(&mut state, &inputs)
            .map_err(AcousticRealizeError::Duct)?;
        self.propagator.sync_nodal(&mut state).map_err(AcousticRealizeError::Duct)?;
        let ledger = &output.ledger;
        if output.port_outputs.len() != inputs.len()
            || output.port_outputs.iter().any(|port|
                !port.conjugate_rate.is_finite() || !port.delta_x.is_finite()
                    || !port.end_conjugate.is_finite())
            || ![
                state.stored_energy_j(), ledger.stored_energy_j,
                ledger.port_work_step_j, ledger.cumulative_port_work_j,
                ledger.cumulative_dissipated_j, ledger.energy_balance_residual_j,
                ledger.relative_energy_drift,
            ].iter().all(|value| value.is_finite())
            || state.pressure().iter().chain(state.volume_flow()).any(|x| !x.is_finite())
        {
            return Err(invalid("DG port candidate left the finite pressure/state/work domain"));
        }
        let inlet_port = output.port_outputs.iter().find(|port| port.id == self.inlet)
            .ok_or_else(|| invalid("DG candidate omitted its original inlet"))?;
        let inlet_work_j = volume_flow_m3_s * inlet_port.delta_x;
        let terminal_frame = if let Some(relation) = relation {
            let outlet_port = output.port_outputs.iter().find(|port| port.id == relation.outlet)
                .ok_or_else(|| invalid("DG candidate omitted its original outlet"))?;
            let h = &relation.response.history_pressure_pa;
            let z = &relation.response.impedance_pa_s_m3;
            let row = relation.outlet_index;
            let pressure_scale = h[row].abs()
                + (z[2 * row + relation.inlet_index] * volume_flow_m3_s).abs()
                + (z[2 * row + row] * outward_flow).abs() + load_pressure[0].abs();
            let port_work_j = outward_flow * outlet_port.delta_x;
            let load_work_j = load_pressure[0] * (outward_flow * self.propagator.dt());
            let port_work_residual_j = port_work_j - load_work_j;
            let work_scale = pressure_scale * outward_flow.abs() * self.propagator.dt();
            if ![
                outward_flow, load_pressure[0], pressure_scale, port_work_j,
                load_work_j, port_work_residual_j, work_scale,
            ].iter().all(|value| value.is_finite())
                || (outlet_port.conjugate_rate - load_pressure[0]).abs()
                    > 1e-10 * pressure_scale.max(1e-300)
                || port_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
            {
                return Err(invalid("DG and rational outlet disagree on same-step pressure or work"));
            }
            Some(DgTerminalFrame {
                port: relation.outlet, outward_flow_m3_s: outward_flow,
                mean_pressure_pa: outlet_port.conjugate_rate,
                port_work_j, load_work_j, port_work_residual_j,
            })
        } else { None };
        let terminal_work = terminal_frame.as_ref().map_or(0.0, |frame| frame.port_work_j);
        let balance = ledger.port_work_step_j - (inlet_work_j - terminal_work);
        let scale = ledger.port_work_step_j.abs() + inlet_work_j.abs() + terminal_work.abs();
        if ![inlet_work_j, balance, scale].iter().all(|value| value.is_finite())
            || balance.abs() > 1e-10 * scale.max(1e-300)
        {
            return Err(invalid("DG original net port work does not equal inlet minus terminal work"));
        }
        Ok(DgPortTrial { owner: self, state, output, terminal, terminal_frame, inlet_work_j })
    }
}

/// A trial borrows the only accepted histories. Dropping it publishes nothing.
pub(crate) struct DgPortTrial<'a> {
    owner: &'a mut DgPortLoad,
    state: DgState,
    output: DgStepOutput,
    terminal: Option<PreparedMimoLoad>,
    terminal_frame: Option<DgTerminalFrame>,
    inlet_work_j: f64,
}

impl DgPortTrial<'_> {
    pub(crate) const fn previous_state(&self) -> &DgState { &self.owner.state }
    pub(crate) const fn state(&self) -> &DgState { &self.state }
    pub(crate) const fn output(&self) -> &DgStepOutput { &self.output }
    pub(crate) const fn inlet_work_j(&self) -> f64 { self.inlet_work_j }
    pub(crate) const fn terminal_frame(&self) -> Option<&DgTerminalFrame> {
        self.terminal_frame.as_ref()
    }

    /// Publish every already validated candidate; all fallible checks precede this.
    pub(crate) fn accept(self) -> DgStepOutput {
        self.owner.state = self.state;
        self.owner.terminal = self.terminal;
        self.owner.last_terminal_frame = self.terminal_frame;
        self.output
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod terminal_tests;
