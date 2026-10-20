//! Simultaneous impedance feedback at explicitly selected network terminals.
//!
//! At each port, `mean(p) = 2a - Zc u = mean(C)x + mean(D)u`. One dense solve
//! couples every declared outlet before any line or rational state advances.
//! Pressure is the exact mean for that held-flow interval, so its product with
//! flow and duration is the rational model's integrated port work. The wave
//! interconnection remains sampled and still requires physical clock refinement.
//! The rational states belong to `MimoZoh`; this module owns only the port
//! interconnection. Signed work delivered to that model is not relabeled as
//! dissipation or reconstructed into a fictitious stored-energy coordinate.

use super::{NetworkNode, NodeFrame, WaveguideError, WaveguideNetwork};
use crate::mimo::{MimoModel, MimoPassivityReport, MimoZoh, check_mimo_passivity};
use crate::model::PoleTerm;
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
    source: MimoTerminalSpec,
    runtime: MimoZoh,
    factor: LuWorkspace,
    impedance: Vec<f64>,
    continuous: MimoPassivityReport,
    sampled: MimoPassivityReport,
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
    fn new(source: MimoTerminalSpec, graph: &WaveguideNetwork) -> Result<Self, WaveguideError> {
        let n = source.model.ports;
        let terms = source.model.poles.len();
        if !(1..=128).contains(&n)
            || source.nodes.len() != n
            || terms > 256
            || source.max_states > 256
            || !(17..=65_537).contains(&source.sampled_passivity_points)
        {
            return Err(invalid(
                "coupled terminals exceed port, pole or sampling limits",
            ));
        }
        let states = source
            .model
            .order()
            .checked_mul(n)
            .ok_or_else(|| invalid("MIMO state size overflow"))?;
        let entries = n * n;
        if states > source.max_states
            || source.model.d.len() != entries
            || source.model.e.len() != entries
            || source.model.residues.len() != terms
            || source.model.residues.iter().any(|r| {
                r.len() != entries || r.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
            || source.model.d.iter().any(|v| !v.is_finite())
            || source.model.e.iter().any(|v| *v != 0.0)
            || !source.model.is_stable()
            || source
                .model
                .poles
                .iter()
                .zip(&source.model.residues)
                .any(|(p, r)| match p {
                    PoleTerm::Real { pole, .. } => {
                        !pole.is_finite() || r.iter().any(|z| z.im != 0.0)
                    }
                    PoleTerm::Pair { pole, .. } => {
                        !pole.re.is_finite() || !pole.im.is_finite() || pole.im <= 0.0
                    }
                })
        {
            return Err(invalid(
                "coupled terminals require a finite stable proper MIMO model with exact shapes",
            ));
        }
        let band = source.band_rad_s;
        if !band.iter().all(|v| v.is_finite())
            || band[0] <= 0.0
            || band[1] <= band[0]
            || band[1] < core::f64::consts::PI / graph.time_step_s
        {
            return Err(invalid(
                "MIMO operator band must cover the physical clock's Nyquist limit",
            ));
        }
        // Overestimate retained matrix/state/record payload before allocating.
        // Finite dense limits above make this arithmetic representable on wasm32.
        let required = entries * (64 * terms + 64)
            + states * 64
            + n * 256
            + (source.sampled_passivity_points + 6400 + terms * 128) * 8;
        if required > source.max_memory_bytes {
            return Err(invalid(
                "coupled terminal payload exceeds its separate memory allowance",
            ));
        }
        let mut impedance = Vec::with_capacity(n);
        for (i, &node) in source.nodes.iter().enumerate() {
            if source.nodes[..i].contains(&node)
                || !matches!(graph.nodes.get(node), Some(NetworkNode::Termination { .. }))
            {
                return Err(invalid(
                    "MIMO ports must be distinct existing memoryless terminals",
                ));
            }
            let range = graph.offsets[node]..graph.offsets[node + 1];
            if range.len() != 1 {
                return Err(invalid("MIMO port must have exactly one adjacent section"));
            }
            let port = graph.ports[range.start];
            impedance.push(graph.lines[port / 2].spec.impedance_pa_s_m3);
        }
        let continuous = check_mimo_passivity(&source.model, (band[0], band[1]))
            .map_err(|_| invalid("MIMO continuous passivity admission failed"))?;
        if !continuous.passive {
            return Err(invalid(
                "MIMO continuous impedance has a passivity violation",
            ));
        }
        let runtime = source
            .model
            .zoh(graph.time_step_s)
            .map_err(|_| invalid("MIMO exact ZOH realization failed"))?;
        let sampled = runtime.mean_sampled_passivity(source.sampled_passivity_points);
        if !sampled.passive {
            return Err(invalid(
                "actual sampled MIMO impedance has a passivity violation",
            ));
        }
        let mut system = runtime.mean_direct_term().to_vec();
        for (i, &z) in impedance.iter().enumerate() {
            system[i * n + i] += z;
        }
        let mut factor =
            LuWorkspace::new(n).map_err(|_| invalid("coupled terminal solve allocation failed"))?;
        factor
            .factor(&system)
            .map_err(|_| invalid("coupled terminal direct solve is singular or nonfinite"))?;
        Ok(Self {
            source,
            runtime,
            factor,
            impedance,
            continuous,
            sampled,
            history: vec![0.0; n],
            rhs: vec![0.0; n],
            flow: vec![0.0; n],
            candidate: vec![MimoPortFrame::default(); n],
            accepted: vec![MimoPortFrame::default(); n],
        })
    }

    pub(super) fn index(&self, node: usize) -> Option<usize> {
        self.source.nodes.iter().position(|&n| n == node)
    }

    pub(super) fn preview(&mut self, arriving: &[f64], dt: f64) -> Result<f64, WaveguideError> {
        let n = self.source.nodes.len();
        self.runtime.mean_history_output_into(&mut self.history);
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
        if self.accepted_steps != 0 || self.mimo.is_some() {
            return Err(invalid(
                "MIMO terminals require an unadvanced network with no existing coupled group",
            ));
        }
        let model = MimoTerminals::new(source, &self)?;
        self.mimo_arriving = vec![0.0; model.source.nodes.len()];
        self.mimo = Some(model);
        Ok(())
    }
    /// Original ordered model, band and work limits actually installed.
    #[must_use]
    pub fn mimo_terminal_spec(&self) -> Option<&MimoTerminalSpec> {
        self.mimo.as_ref().map(|m| &m.source)
    }
    /// Continuous and actual sampled passivity reports, without conflating them.
    #[must_use]
    pub fn mimo_passivity(&self) -> Option<(&MimoPassivityReport, &MimoPassivityReport)> {
        self.mimo.as_ref().map(|m| (&m.continuous, &m.sampled))
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
        for (i, &node) in model.source.nodes.iter().enumerate() {
            let port = self.ports[self.offsets[node]];
            self.mimo_arriving[i] = self.arriving(port);
        }
        let work = self
            .mimo
            .as_mut()
            .expect("same admitted group")
            .preview(&self.mimo_arriving, self.time_step_s)?;
        let model = self.mimo.as_ref().expect("same admitted group");
        for (i, &node) in model.source.nodes.iter().enumerate() {
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
