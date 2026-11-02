//! Simultaneous impedance feedback at explicitly selected network terminals.
//!
//! At each port, `mean(p) = 2a - Zc u = mean(C)x + mean(D)u`. One dense solve
//! couples every declared outlet before any line or rational state advances.
//! Pressure is the exact mean for that held-flow interval, so its product with
//! flow and duration is the rational model's integrated port work. The wave
//! interconnection remains sampled and still requires physical clock refinement.
//! Pole and propagation histories belong to the original admitted runtime;
//! this module owns only the port
//! interconnection. Signed work delivered to that model is not relabeled as
//! dissipation or reconstructed into a fictitious stored-energy coordinate.

use super::{NetworkNode, NodeFrame, WaveguideError, WaveguideNetwork};
use crate::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use crate::mimo::{MimoModel, MimoPassivityReport, MimoZoh};
use crate::mimo::delayed::{DelayedMimoLoadLimits, PreparedDelayedMimoLoad};
use crate::mimo::DelayedMimoModel;
use fs_la::prepared_lu::LuWorkspace;

/// An explicitly ordered, proper rational terminal impedance in SI units.
#[derive(Clone, Debug)]
pub struct MimoTerminalSpec {
    /// Existing degree-one memoryless nodes, in model row/column order.
    pub nodes: Vec<usize>,
    /// `p[Pa] = Z[Pa s/m^3] u[m^3/s]`, on the `s = +i omega` axis.
    /// No per-entry propagation delays or nonzero improper `s E` terms.
    pub model: MimoModel,
    /// Declared operator band [rad/s]. It must cover the source Nyquist limit.
    pub band_rad_s: [f64; 2],
    /// Actual mean-output ZOH check, from DC through Nyquist; 17..=65537.
    /// More samples strengthen a diagnostic, never turn it into a certificate.
    pub sampled_passivity_points: usize,
    /// Maximum real state dimension, bounded by 256 in this dense adapter.
    pub max_states: usize,
    /// Additional retained payload allowance, separate from original wave memory.
    /// Excludes input vectors, allocator overhead and cold eigensolver scratch.
    pub max_memory_bytes: usize,
}

/// An ordered terminal impedance retaining every physical propagation delay.
/// The source band, physical clock and separate pole/history/memory limits are
/// admitted by the delayed owner. Both response checks remain `GridOnly`.
#[derive(Clone, Debug)]
pub struct DelayedMimoTerminalSpec {
    /// Distinct degree-one memoryless nodes, in original matrix order.
    pub nodes: Vec<usize>,
    /// Proper impedance on `s = +i omega`, with per-entry delays in seconds.
    pub model: DelayedMimoModel,
    /// Runtime limits; its memory allowance also covers this terminal adapter.
    pub limits: DelayedMimoLoadLimits,
}

enum TerminalSource {
    Rational(MimoTerminalSpec),
    Delayed(DelayedMimoTerminalSpec),
}
impl TerminalSource {
    fn nodes(&self) -> &[usize] {
        match self { Self::Rational(s) => &s.nodes, Self::Delayed(s) => &s.nodes }
    }
    fn ports(&self) -> usize {
        match self { Self::Rational(s) => s.model.ports, Self::Delayed(s) => s.model.model.ports }
    }
}

enum TerminalRuntime {
    Rational {
        runtime: MimoZoh,
        continuous: MimoPassivityReport,
        sampled: MimoPassivityReport,
    },
    Delayed(PreparedDelayedMimoLoad),
}
impl TerminalRuntime {
    fn mean_direct_term(&self) -> &[f64] {
        match self { Self::Rational { runtime, .. } => runtime.mean_direct_term(), Self::Delayed(r) => r.mean_direct_term() }
    }
    fn mean_history_output_into(&self, output: &mut [f64]) -> Result<(), WaveguideError> {
        match self {
            Self::Rational { runtime, .. } => { runtime.mean_history_output_into(output); Ok(()) }
            Self::Delayed(r) => r.mean_history_output_into(output).map_err(|e| invalid(e.0)),
        }
    }
    fn preview_state(&mut self, input: &[f64]) -> bool {
        match self { Self::Rational { runtime, .. } => runtime.preview_state(input), Self::Delayed(r) => r.preview_state(input) }
    }
    fn accept_state(&mut self, input: &[f64]) {
        match self { Self::Rational { runtime, .. } => runtime.accept_state(input), Self::Delayed(r) => r.accept_state(input) }
    }
    fn pole_states(&self) -> &[fs_math::c64::C64] {
        match self { Self::Rational { runtime, .. } => runtime.pole_states(), Self::Delayed(r) => r.pole_states() }
    }
    fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        match self {
            Self::Rational { continuous, sampled, .. } => (continuous, sampled),
            Self::Delayed(runtime) => runtime.reports(),
        }
    }
}

/// Accepted same-sample effort, flow and signed work at one coupled port.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MimoPortFrame {
    /// Interval-mean exterior pressure returned to this aperture [Pa].
    pub pressure_pa: f64,
    /// Signed volume flow from the network into the exterior [m^3/s].
    pub flow_m3_s: f64,
    /// Exact held-flow interval work `mean(pressure) * flow * dt` [J].
    /// May be negative as the model returns energy.
    pub work_j: f64,
}

pub(super) struct MimoTerminals {
    source: TerminalSource,
    runtime: TerminalRuntime,
    factor: LuWorkspace,
    impedance: Vec<f64>,
    history: Vec<f64>,
    rhs: Vec<f64>,
    flow: Vec<f64>,
    candidate: Vec<MimoPortFrame>,
    accepted: Vec<MimoPortFrame>,
}

fn invalid(what: &'static str) -> WaveguideError {
    WaveguideError(what)
}

impl MimoTerminals {
    fn new(source: TerminalSource, graph: &WaveguideNetwork) -> Result<Self, WaveguideError> {
        let n = source.ports();
        if !(1..=128).contains(&n) || source.nodes().len() != n {
            return Err(invalid("coupled terminals require one existing node per model port"));
        }
        let mut impedance = Vec::with_capacity(n);
        for (i, &node) in source.nodes().iter().enumerate() {
            if source.nodes()[..i].contains(&node)
                || !matches!(graph.nodes.get(node), Some(NetworkNode::Termination { .. }))
            {
                return Err(invalid("MIMO ports must be distinct existing memoryless terminals"));
            }
            let range = graph.offsets[node]..graph.offsets[node + 1];
            if range.len() != 1 {
                return Err(invalid("MIMO port must have exactly one adjacent section"));
            }
            let port = graph.ports[range.start];
            impedance.push(graph.lines[port / 2].spec.impedance_pa_s_m3);
        }
        let runtime = match &source {
            TerminalSource::Rational(source) => {
                let limits = MimoLoadLimits {
                    band_rad_s: source.band_rad_s,
                    sampled_passivity_points: source.sampled_passivity_points,
                    max_states: source.max_states,
                    max_memory_bytes: source.max_memory_bytes,
                };
                let (runtime, continuous, sampled) =
                    PreparedMimoLoad::new(&source.model, graph.time_step_s, limits)
                        .map_err(|error| invalid(error.0))?.into_parts();
                TerminalRuntime::Rational { runtime, continuous, sampled }
            }
            TerminalSource::Delayed(source) => {
                let terms = source.model.model.poles.len();
                if terms > 256 { return Err(invalid("MIMO load exceeds port, pole or sampling limits")); }
                // The delayed owner retains its source. Charge the original
                // installed specification and this adapter's dense solve,
                // vectors and port frames before granting its remaining budget.
                let retained = n * n * (16 * terms + 24 + 32)
                    + terms * core::mem::size_of::<crate::model::PoleTerm>()
                    + n * (80 + 2 * core::mem::size_of::<MimoPortFrame>());
                let mut limits = source.limits;
                limits.source.max_memory_bytes = limits.source.max_memory_bytes
                    .checked_sub(retained)
                    .ok_or_else(|| invalid("delayed MIMO terminal adapter exceeds its memory allowance"))?;
                let runtime = PreparedDelayedMimoLoad::new(&source.model, graph.time_step_s, limits)
                    .map_err(|error| invalid(error.0))?;
                TerminalRuntime::Delayed(runtime)
            }
        };
        let mut system = runtime.mean_direct_term().to_vec();
        for (i, &z) in impedance.iter().enumerate() {
            system[i * n + i] += z;
        }
        let mut factor =
            LuWorkspace::new(n).map_err(|_| invalid("coupled terminal solve allocation failed"))?;
        factor.factor(&system)
            .map_err(|_| invalid("coupled terminal direct solve is singular or nonfinite"))?;
        Ok(Self {
            source, runtime, factor, impedance,
            history: vec![0.0; n], rhs: vec![0.0; n], flow: vec![0.0; n],
            candidate: vec![MimoPortFrame::default(); n],
            accepted: vec![MimoPortFrame::default(); n],
        })
    }

    pub(super) fn index(&self, node: usize) -> Option<usize> {
        self.source.nodes().iter().position(|&n| n == node)
    }

    pub(super) fn preview(&mut self, arriving: &[f64], dt: f64) -> Result<f64, WaveguideError> {
        let n = self.source.nodes().len();
        self.runtime.mean_history_output_into(&mut self.history)?;
        for i in 0..n {
            self.rhs[i] = 2.0 * arriving[i] - self.history[i];
        }
        self.factor
            .solve_factored_into(&self.rhs, &mut self.flow)
            .map_err(|_| invalid("coupled terminal solve left the finite set"))?;
        let mut work = 0.0;
        for i in 0..n {
            let pressure = self.history[i]
                + self.runtime.mean_direct_term()[i * n..(i + 1) * n]
                    .iter()
                    .zip(&self.flow)
                    .map(|(d, u)| d * u)
                    .sum::<f64>();
            let expected = 2.0 * arriving[i] - self.impedance[i] * self.flow[i];
            let scale = pressure.abs()
                + self.history[i].abs()
                + 2.0 * arriving[i].abs()
                + (self.impedance[i] * self.flow[i]).abs();
            let port_work = pressure * self.flow[i] * dt;
            if ![pressure, expected, self.flow[i], port_work, scale]
                .iter()
                .all(|v| v.is_finite())
                || (pressure - expected).abs() > 1024.0 * f64::EPSILON * scale
            {
                return Err(invalid(
                    "MIMO pressure/flow solve is nonfinite or exceeds its roundoff budget",
                ));
            }
            self.candidate[i] = MimoPortFrame {
                pressure_pa: pressure,
                flow_m3_s: self.flow[i],
                work_j: port_work,
            };
            work += port_work;
        }
        if !work.is_finite() || !self.runtime.preview_state(&self.flow) {
            return Err(invalid(
                "MIMO candidate state or total work left the finite set",
            ));
        }
        Ok(work)
    }
    pub(super) fn frame(&self, i: usize) -> MimoPortFrame {
        self.candidate[i]
    }
    pub(super) fn accept(&mut self) {
        self.runtime.accept_state(&self.flow);
        self.accepted.copy_from_slice(&self.candidate);
    }
}

impl WaveguideNetwork {
    /// Install one coupled terminal group before any accepted sample. Existing
    /// reactive boundaries, interior loads, geometry and wave buffers are kept.
    /// A second group or any later attachment refuses; no history is discarded.
    ///
    /// The continuous certificate is retained; actual mean-output ZOH response
    /// is checked independently on the declared full-Nyquist grid. Its `GridOnly` result
    /// does not establish coupled discrete passivity or rational stored energy.
    /// The network's balance therefore accounts for signed external work at the
    /// MIMO boundary separately from its own storage and irreversible losses.
    pub fn with_mimo_terminals(mut self, source: MimoTerminalSpec) -> Result<Self, WaveguideError> {
        self.install_mimo_terminals(source)?;
        Ok(self)
    }
    /// Cold in-place form of `with_mimo_terminals`. Refusal leaves the original
    /// network's accepted physics and terminal ownership unchanged.
    pub fn install_mimo_terminals(
        &mut self,
        source: MimoTerminalSpec,
    ) -> Result<(), WaveguideError> {
        self.install_terminal_group(TerminalSource::Rational(source))
    }
    /// Install an exact held-input delayed impedance before any accepted sample.
    /// Fractional delays split the physical interval analytically; no delay
    /// rounding, interpolation filter or additional feedback sample is inserted.
    /// Pressure, flow and signed work use the same simultaneous terminal solve.
    /// The restored continuous and actual discrete checks remain `GridOnly`;
    /// they do not establish a delayed storage functional or closed-loop stability.
    pub fn with_delayed_mimo_terminals(
        mut self, source: DelayedMimoTerminalSpec,
    ) -> Result<Self, WaveguideError> {
        self.install_delayed_mimo_terminals(source)?;
        Ok(self)
    }
    /// Failure-atomic cold attachment of one explicitly delayed terminal group.
    pub fn install_delayed_mimo_terminals(
        &mut self, source: DelayedMimoTerminalSpec,
    ) -> Result<(), WaveguideError> {
        self.install_terminal_group(TerminalSource::Delayed(source))
    }
    fn install_terminal_group(&mut self, source: TerminalSource) -> Result<(), WaveguideError> {
        if self.accepted_steps != 0 || self.mimo.is_some() {
            return Err(invalid(
                "MIMO terminals require an unadvanced network with no existing coupled group",
            ));
        }
        let model = MimoTerminals::new(source, &self)?;
        self.mimo_arriving = vec![0.0; model.source.nodes().len()];
        self.mimo = Some(model);
        Ok(())
    }
    /// Original undelayed model, band and limits; absent for a delayed group.
    #[must_use]
    pub fn mimo_terminal_spec(&self) -> Option<&MimoTerminalSpec> {
        self.mimo.as_ref().and_then(|m| match &m.source {
            TerminalSource::Rational(source) => Some(source), TerminalSource::Delayed(_) => None,
        })
    }
    /// Original delayed model, node order and complete caller limits.
    #[must_use]
    pub fn delayed_mimo_terminal_spec(&self) -> Option<&DelayedMimoTerminalSpec> {
        self.mimo.as_ref().and_then(|m| match &m.source {
            TerminalSource::Delayed(source) => Some(source), TerminalSource::Rational(_) => None,
        })
    }
    /// Read-only accepted delayed poles, input rings and reports. Its memory
    /// allowance is the remainder after the terminal adapter's retained payload.
    #[must_use]
    pub fn delayed_mimo_load(&self) -> Option<&PreparedDelayedMimoLoad> {
        self.mimo.as_ref().and_then(|m| match &m.runtime {
            TerminalRuntime::Delayed(runtime) => Some(runtime), TerminalRuntime::Rational { .. } => None,
        })
    }
    /// Continuous and actual sampled passivity reports, without conflating them.
    #[must_use]
    pub fn mimo_passivity(&self) -> Option<(&MimoPassivityReport, &MimoPassivityReport)> {
        self.mimo.as_ref().map(|m| m.runtime.reports())
    }
    /// Last accepted coupled port observations, in original model order.
    #[must_use]
    pub fn mimo_port_frames(&self) -> Option<&[MimoPortFrame]> {
        self.mimo.as_ref().map(|m| m.accepted.as_slice())
    }
    /// Accepted original ZOH pole coordinates; failed samples cannot replace them.
    #[must_use]
    pub fn mimo_pole_states(&self) -> Option<&[fs_math::c64::C64]> {
        self.mimo.as_ref().map(|m| m.runtime.pole_states())
    }

    pub(super) fn preview_mimo(&mut self) -> Result<f64, WaveguideError> {
        let Some(model) = self.mimo.as_ref() else {
            return Ok(0.0);
        };
        for (i, &node) in model.source.nodes().iter().enumerate() {
            let port = self.ports[self.offsets[node]];
            self.mimo_arriving[i] = self.arriving(port);
        }
        let work = self
            .mimo
            .as_mut()
            .expect("same admitted group")
            .preview(&self.mimo_arriving, self.time_step_s)?;
        let model = self.mimo.as_ref().expect("same admitted group");
        for (i, &node) in model.source.nodes().iter().enumerate() {
            let port = self.ports[self.offsets[node]];
            let frame = model.frame(i);
            self.departing[port] = frame.pressure_pa - self.mimo_arriving[i];
            if !self.departing[port].is_finite() {
                return Err(invalid("MIMO returning wave is not finite"));
            }
            self.candidate[node] = NodeFrame {
                pressure_pa: frame.pressure_pa,
                net_flow_into_node_m3_s: frame.flow_m3_s,
                load_pressure_pa: frame.pressure_pa,
                load_flow_m3_s: frame.flow_m3_s,
                ..NodeFrame::default()
            };
        }
        Ok(work)
    }
}

#[cfg(test)]
mod tests;
