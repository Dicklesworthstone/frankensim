//! A physical tube graph under the existing jet-drive island.
//!
//! No second jet equation or propagation integrator. Original node addresses
//! survive viscothermal subdivision, so loads, branches and outlet observations
//! refer to the same geometry. Gas fields are frozen, common-pressure states;
//! mixing, mean-flow convection and time-varying materials are not modeled.

use super::{FlueError, FlueGeometry, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{SectionRealization, TubeNetworkSpec};
use crate::bernoulli_aperture::viscothermal::{
    ViscothermalSection, ViscothermalSectionSpec, with_regional_viscothermal_sections,
};
use crate::wind_line::WindLineBank;
use fs_aeroac::jetcard::JetCard;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::waveguide::network::{NetworkFrame, NodeFrame, WaveguideNetwork};

/// Original geometry, state and explicit approximation/work choices.
#[derive(Clone, Debug)]
pub struct FlueNetworkSpec {
    /// Original cylinders, branches and passive load declarations.
    pub graph: TubeNetworkSpec,
    /// One complete material-gas state per ORIGINAL section.
    pub gases: Vec<GasState>,
    /// Selected distributed viscous/thermal wall laws; empty means no such loss.
    /// Existing wide-tube applicability and sampled-error gates are unchanged.
    pub viscothermal: Vec<ViscothermalSectionSpec>,
    /// Lowest strictly positive jet speed [m/s] admitted during performance.
    /// Zero blowing is independently admitted and preserves stored vibration.
    pub minimum_jet_speed_m_s: f64,
    /// Explicit incompressible-jet Mach ceiling, in (0, 0.3].
    pub maximum_jet_mach: f64,
    /// Maximum retained jet convection samples, including interpolation slots.
    /// Separate from graph payload; bounded above by 1,048,576 samples.
    pub max_jet_history_samples: usize,
}

/// The existing wave/load owner and its exact original material/geometry inputs.
/// The jet card's estimate-only authority is not changed by this construction.
pub struct FlueNetwork {
    source: FlueNetworkSpec,
    graph: TubeNetworkSpec,
    gases: Vec<GasState>,
    represented: Vec<SectionRealization>,
    losses: Vec<ViscothermalSection>,
    physics: WaveguideNetwork,
    last: Option<NetworkFrame>,
}

pub(super) enum FlueResonator {
    Tmm(WindLineBank),
    Network(Box<FlueNetwork>),
}

fn invalid(what: &'static str) -> FlueError { FlueError::Invalid { what } }
fn physical(e: AcousticRealizeError) -> FlueError { FlueError::Physics(e) }
fn wave(e: fs_vfit::waveguide::WaveguideError) -> FlueError {
    physical(AcousticRealizeError::Nonlinear(e.to_string()))
}
fn poll(gate: &CancelGate) -> Result<(), FlueError> {
    if gate.is_requested() { Err(physical(AcousticRealizeError::Cancelled)) } else { Ok(()) }
}

impl FlueNetwork {
    fn new(source: FlueNetworkSpec, dt: f64, zc: f64, gate: &CancelGate) -> Result<Self, FlueError> {
        poll(gate)?;
        let (graph, gases, losses) = with_regional_viscothermal_sections(
            source.graph.clone(), source.gases.clone(), dt, &source.viscothermal, gate,
        ).map_err(physical)?;
        // A caller cannot provide a second, inconsistent source cross-section.
        let actual = graph.inlet_impedance_with_gases(&gases).map_err(physical)?;
        if actual.to_bits() != zc.to_bits() {
            return Err(invalid("flue mouth bore area must match the physical inlet section"));
        }
        let (physics, represented) = graph.realize_waveguide(&gases, dt, gate).map_err(physical)?;
        poll(gate)?;
        Ok(Self { source, graph, gases, represented, losses, physics, last: None })
    }
    /// Original data before any loss-cell subdivision.
    #[must_use]
    pub const fn source(&self) -> &FlueNetworkSpec { &self.source }
    /// Actual lowered graph. Original node indices remain valid.
    #[must_use]
    pub const fn graph(&self) -> &TubeNetworkSpec { &self.graph }
    /// Complete ordered states actually used for propagation and observation.
    #[must_use]
    pub fn section_gases(&self) -> &[GasState] { &self.gases }
    /// Physical lengths and impedances admitted by the shared geometry owner.
    #[must_use]
    pub fn represented_sections(&self) -> &[SectionRealization] { &self.represented }
    /// Source-derived viscothermal spectra and sampled discretization errors.
    #[must_use]
    pub fn viscothermal_sections(&self) -> &[ViscothermalSection] { &self.losses }
    /// Last accepted node pressure/flow; initially zero.
    #[must_use]
    pub fn node_frame(&self, node: usize) -> Option<&NodeFrame> { self.physics.node_frame(node) }
    /// Last accepted whole-network work and storage, absent before stepping.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&NetworkFrame> { self.last.as_ref() }
    /// Actual wave plus reactive-load storage, not an integrated work estimate.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 { self.physics.stored_energy_j() }
}

impl FlueResonator {
    pub(super) fn incoming(&self) -> f64 {
        match self { Self::Tmm(line) => line.incoming(), Self::Network(n) => n.physics.incoming_pressure_pa() }
    }
    pub(super) fn push(&mut self, outgoing: f64) -> Result<f64, FlueError> {
        match self {
            Self::Tmm(line) => line.push(outgoing).map_err(FlueError::Line),
            Self::Network(n) => {
                let frame = n.physics.step(outgoing).map_err(wave)?;
                n.last = Some(frame);
                Ok(n.physics.incoming_pressure_pa())
            }
        }
    }
    pub(super) fn switch_fingering(&mut self, index: usize, fade: usize) -> Result<(), FlueError> {
        match self {
            Self::Tmm(line) => line.switch_fingering(index, fade).map_err(FlueError::Line),
            Self::Network(_) => Err(invalid("a physical flue graph has no TMM fingering bank; use explicit terminal controls")),
        }
    }
    pub(super) fn validate_pressure(&self, p: f64, rho: f64) -> Result<(), FlueError> {
        if let Self::Network(n) = self {
            let source = &n.source;
            let speed = (2.0 * p / rho).sqrt();
            let inlet = source.graph.inlet_section_index().map_err(physical)?;
            if !speed.is_finite() || (p > 0.0 && speed < source.minimum_jet_speed_m_s)
                || speed > source.maximum_jet_mach * source.gases[inlet].sound_speed {
                return Err(invalid("blowing pressure exceeds the declared convection-history or low-Mach jet domain"));
            }
        }
        Ok(())
    }
}

impl FlueVoice {
    /// Compose the ORIGINAL jet-drive loop with physical regional tube sections.
    ///
    /// Geometry sets delays/impedances; selected viscothermal laws and supplied
    /// passive loads supply their own stored states and work. The same inlet gas
    /// drives Bernoulli speed, jet Reynolds number, convection and mouth loading.
    /// No fake reed, oscillator frequency, material-name switch or extra sample
    /// of junction delay is introduced. `step_block` still returns internal p+;
    /// exterior pressure requires an explicitly placed flow observer.
    ///
    /// This preserves the existing tonal card, authored receptivity gain and
    /// seed: it does NOT claim turbulent flute noise or empirical instrument
    /// accuracy. Ideal junctions omit branch matching/end corrections unless
    /// these are supplied as explicit physical loads. No mean gas flow is solved.
    #[allow(clippy::too_many_arguments)]
    pub fn from_network(
        card: Option<&JetCard>, displacement_gain: f64, seed_amplitude_rel: f64,
        geometry: FlueGeometry, source: FlueNetworkSpec, sample_rate_hz: u32,
        voice_seed: u64, voice_index: u32, gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        poll(gate)?;
        // Preserve the no-card refusal before touching expensive reduction.
        let card = card.ok_or(FlueError::NoCard)?;
        source.graph.validate_section_gases(&source.gases).map_err(physical)?;
        if !source.minimum_jet_speed_m_s.is_finite() || source.minimum_jet_speed_m_s <= 0.0
            || !source.maximum_jet_mach.is_finite() || !(0.0 < source.maximum_jet_mach && source.maximum_jet_mach <= 0.3)
            || !(4..=1_048_576).contains(&source.max_jet_history_samples) {
            return Err(invalid("flue network needs explicit positive jet speed, low-Mach and history limits"));
        }
        let gas = source.gases[source.graph.inlet_section_index().map_err(physical)?];
        if source.minimum_jet_speed_m_s > source.maximum_jet_mach * gas.sound_speed {
            return Err(invalid("flue jet speed interval is empty"));
        }
        let minimum = source.minimum_jet_speed_m_s;
        let capacity = source.max_jet_history_samples;
        let voice = Self::with_resonator(Some(card), displacement_gain, seed_amplitude_rel,
            geometry, &gas, sample_rate_hz, voice_seed, voice_index, minimum, capacity,
            |zc| FlueNetwork::new(source, 1.0 / f64::from(sample_rate_hz), zc, gate)
                .map(|n| FlueResonator::Network(Box::new(n))))?;
        poll(gate)?;
        Ok(voice)
    }
    /// Read the actual regional tube state; absent for the legacy TMM image.
    #[must_use]
    pub fn network(&self) -> Option<&FlueNetwork> {
        match &self.line { FlueResonator::Network(n) => Some(n), FlueResonator::Tmm(_) => None }
    }
    /// Change an ideal memoryless terminal without resetting any traveling wave.
    /// This is not moving-pad contact, displacement work or a radiation model.
    /// Reactive loads refuse rather than discard their stored energy.
    pub fn set_network_terminal_reflection(&mut self, node: usize, reflection: f64) -> Result<(), FlueError> {
        if self.poisoned { return Err(invalid("flue voice is poisoned by a previous failed sample")); }
        match &mut self.line {
            FlueResonator::Network(n) => {
                n.physics.set_terminal_reflection(node, reflection).map_err(wave)?;
                n.graph.nodes[node] = fs_vfit::waveguide::network::NetworkNode::Termination { reflection };
                Ok(())
            }
            FlueResonator::Tmm(_) => Err(invalid("terminal controls require the physical graph image")),
        }
    }
    /// Validate a held physical pressure before compiling a source schedule.
    pub fn validate_blowing_pressure(&self, pressure_pa: f64) -> Result<(), FlueError> {
        if self.poisoned { return Err(invalid("flue voice is poisoned by a previous failed sample")); }
        if !pressure_pa.is_finite() || pressure_pa < 0.0 {
            return Err(invalid("blowing pressure must be finite and non-negative"));
        }
        self.line.validate_pressure(pressure_pa, self.rho)
    }
}

#[cfg(test)]
mod tests;
