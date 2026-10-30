//! Authored pressure phrases on retained coupled valve/tube mechanics.
//!
//! This adapter owns the existing system, not a replacement voice. Plate shape,
//! spatial closure, material memory, wall/load state and traveling waves remain
//! in their physical owners. Pressure and explicit mechanical forces are scheduled. The output
//! is either a named internal pressure or an explicitly selected one-way baffled
//! outlet receiver. The latter observes actual terminal FLOW, not bore pressure.

/// Supplied mesh, closure, material memory and canonical pressure sources.
pub mod file;
/// Physical force programs on the original plate, independent of blowing pressure.
pub mod force;
/// Explicit pressures and forces on the retained moving graph boundaries.
pub mod boundaries;
/// Scheduled heat into the original finite thermochemical DG inventory.
pub mod finite_gas;
/// Conservative finite reservoir exchange on the original source clock.
pub mod finite_transport;
/// Scheduled finite-body wall heating through the original coupled DG owner.
pub mod thermal_wall;

use super::dynamic::DynamicAperture;
use super::dg::ApertureDg;
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use super::network::ApertureNetwork;
use super::network::boundaries::BoundaryNetworkFrame;
use super::tube::{ApertureTube, TubeDrive};
use crate::pcm_wav::observation::{ObservationError, ObservationPolicy, PressureRenderer, validate_reaction_medium};
use super::radiation::BaffledRadiationLoad;
use crate::render::{ControlDelta, RenderError};
use crate::render::schedule::{
    GestureCompileError, PressureGestureBinding, ScheduledControl, compile_pressure_gestures,
};
use fs_scenario::gesture::GestureSchedule;
use fs_material::gas::GasState;

/// Existing reciprocal physical systems; both retain their accepted state.
pub enum CoupledAperture {
    /// One geometry-bound uniform tube and its termination.
    Tube(ApertureTube),
    /// Explicit connected sections, junctions and any already-attached loads.
    Network(ApertureNetwork),
    /// Exact DG acoustic/wall memory coupled through the original mean-flow valve.
    Dg(ApertureDg),
}
impl CoupledAperture {
    /// Immutable access to the actual valve, including plate/contact/material data.
    #[must_use]
    pub fn aperture(&self) -> &DynamicAperture {
        match self { Self::Tube(m) => m.aperture(), Self::Network(m) => m.aperture(),
            Self::Dg(m) => m.aperture() }
    }
}

/// Which physical pressure is sampled once per complete mechanical step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ApertureObservation {
    /// Junction pressure used in the implicit pressure/structure solve.
    Inlet,
    /// Actual terminal pressure after propagation through a uniform tube.
    TubeTerminal,
    /// Actual accepted node pressure in a connected tube network.
    NetworkNode(usize),
    /// Exterior receiver of a tube outlet or every bound DG terminal flow.
    /// The circular openings lie in an infinite rigid baffle. A legacy tube keeps
    /// its authored reflection; a DG outlet requires its computed baffled load
    /// and uses that load's original geometry and exterior gas. Observation
    /// adds no second reaction. DG uses the existing same-clock differentiating
    /// image; the legacy tube's finite-difference image remains unchanged.
    TubeBaffled(CircularOutletReceiver),
    /// Exterior observation of flow into a degree-one passive network terminal.
    /// The adjacent section supplies the actual outlet radius. Its existing load
    /// determines reaction; observation does not add a second load or assume that
    /// an arbitrary supplied impedance represents this baffled mouth.
    NetworkBaffled {
        /// Degree-one passive terminal, not a junction, shunt, series node or inlet.
        node: usize,
        /// Position/quadrature relative to that outlet's outward-facing baffle.
        receiver: CircularOutletReceiver,
    },
}

/// Finite window and bounded offline pressure compilation.
#[derive(Clone, Copy, Debug)]
pub struct AperturePerformanceConfig {
    /// Must exactly match dt = 1 / sample_rate_hz in the supplied physical system.
    pub sample_rate_hz: u32,
    /// Positive half-open source window [0, samples); must fit its step budget.
    pub samples: u64,
    /// Largest source callback; independent of control timing.
    pub max_block: usize,
    /// Bounds canonical pressure compilation, and separately the combined cold
    /// boundary-control admission and source-force projection work.
    pub max_compile_work: u64,
    /// Bounds sampled pressure ticks AND all stored pressure, source-force and
    /// boundary assignments (including explicit initial boundary drives).
    pub max_controls: usize,
}

/// Pressure and optional mechanical-force playback retaining one physical clock.
///
/// Starts with an unadvanced system, which may have nonzero physical initial
/// motion/contact/material energy. The supplied track owns the pressure history;
/// no fixture attack envelope, reset at release, note-to-pressure calibration or
/// prescribed body-flow source is added. Continue by retaining this whole object.
/// Scheduling itself allocates nothing inside callbacks; physical-owner
/// allocations are unchanged. A failed physical callback poisons this adapter,
/// as required by PressureRenderer, even though the owner rolls back its failed
/// individual step. Completed earlier samples remain inspectable.
pub struct AperturePerformance {
    system: CoupledAperture,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    dg_terminal_flows: Vec<f64>,
    exterior_gas: Option<GasState>,
    config: AperturePerformanceConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    poisoned: bool,
    forces: Option<force::ForceSchedule>,
    boundaries: Option<boundaries::BoundarySchedule>,
    boundary_frame: Option<BoundaryNetworkFrame>,
    mechanical_work_j: f64,
    finite_gas: Option<finite_gas::FiniteGasPlayback>,
    finite_transport: Option<finite_transport::FiniteGasTransportPlayback>,
    finite_wall: Option<thermal_wall::FiniteWallPlayback>,
    observation_policy: ObservationPolicy,
}
fn receiver_error(error: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error))
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }

impl AperturePerformance {
    /// Bind acceptance to an actual installed load, including its medium, radius,
    /// clock and resolved band. A numeric reflection is diagnostic-only.
    /// The load is checked against the exclusively owned network before return.
    pub fn new_with_observation_policy(
        system: CoupledAperture,
        observation: ApertureObservation,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
        exterior: Option<GasState>,
        reaction: Option<BaffledRadiationLoad>,
        policy: ObservationPolicy,
    ) -> Result<Self, GestureCompileError> {
        let mut source = Self::new_in_medium(system, observation, schedule, config, exterior)?;
        if policy == ObservationPolicy::Acceptance || reaction.is_some() {
            let load = reaction.ok_or_else(|| {
                GestureCompileError::Render(receiver_error(
                    ObservationError::UnmatchedReactionLoad.to_string(),
                ))
            })?;
            let medium = source.baffled_medium().ok_or_else(|| {
                GestureCompileError::Render(receiver_error(
                    ObservationError::UnmatchedReactionLoad.to_string(),
                ))
            })?;
            validate_reaction_medium(
                RayleighMedium {
                    density: load.density_kg_m3(),
                    sound_speed: load.sound_speed_m_s(),
                },
                medium,
            )
            .map_err(|e| GestureCompileError::Render(receiver_error(e.to_string())))?;
            let matched = match (&source.system, observation) {
                (
                    CoupledAperture::Network(network),
                    ApertureObservation::NetworkBaffled { node, receiver },
                ) => {
                    let mut adjacent = network
                        .spec()
                        .sections
                        .iter()
                        .filter(|s| s.nodes.contains(&node));
                    network.spec().nodes.get(node) == Some(&load.termination())
                        && adjacent
                            .next()
                            .is_some_and(|s| s.radius_m == load.radius_m())
                        && adjacent.next().is_none()
                        && receiver.maximum_frequency_hz <= load.maximum_frequency_hz()
                        && load.time_step_s().to_bits()
                            == (1.0 / f64::from(config.sample_rate_hz)).to_bits()
                }
                _ => false,
            };
            if !matched {
                return Err(GestureCompileError::Render(receiver_error(
                    ObservationError::ReactionBinding {
                        what: "terminal, radius, clock or band differs from the installed load",
                    }
                    .to_string(),
                )));
            }
        }
        source.observation_policy = policy;
        Ok(source)
    }

    /// The observation's admitted use, independent of mechanics or recording fidelity.
    pub const fn observation_policy(&self) -> ObservationPolicy {
        self.observation_policy
    }
    /// Compile and admit the complete pressure phrase without advancing physics.
    /// No source command may begin after the last observed control tick. Ramps
    /// extending beyond the window retain only their actually observed portion.
    ///
    /// # Errors
    /// Clock/horizon/observer mismatch, unsupported or multiple tracks, exhausted
    /// compile/storage budgets, or a source gesture refusal. No runtime is returned
    /// with an incomplete schedule or a substituted pressure observation.
    pub fn new(system: CoupledAperture, observation: ApertureObservation,
        schedule: GestureSchedule, config: AperturePerformanceConfig)
        -> Result<Self, GestureCompileError>
    {
        Self::new_in_medium(system, observation, schedule, config, None)
    }

    /// Observe an unadvanced outlet in an explicitly resolved exterior gas.
    /// The geometry and outward volume flow remain the selected terminal's.
    /// A network selects its receiver medium; a DG source must already bind the
    /// same exterior for its computed radiation reaction. This does not replace
    /// an installed load or change any internal material or physical history.
    ///
    /// A network must retain its complete section gases so common static
    /// pressure can be checked. A DG receiver is already bound to its computed
    /// load's exterior; an explicit exterior must equal that original state.
    /// Legacy scalar-only tubes cannot supply this condition and refuse; use an
    /// explicitly gas-bound one-section network. The original local-gas paths
    /// retain their existing receiver arithmetic.
    ///
    /// # Errors
    /// Original clock/schedule/receiver refusal, incoherent exterior, unequal
    /// pressure, missing section gas binding or a non-exterior observation.
    pub fn new_with_exterior_gas(system: CoupledAperture, observation: ApertureObservation,
        schedule: GestureSchedule, config: AperturePerformanceConfig, exterior: GasState)
        -> Result<Self, GestureCompileError>
    {
        Self::new_in_medium(system, observation, schedule, config, Some(exterior))
    }

    fn new_in_medium(
        system: CoupledAperture,
        observation: ApertureObservation,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
        exterior_gas: Option<GasState>,
    ) -> Result<Self, GestureCompileError> {
        let a = system.aperture();
        if config.sample_rate_hz == 0 || config.samples == 0 || config.max_block == 0
            || a.accepted_steps() != 0 || config.samples > a.spec().max_steps
            || a.spec().time_step_s.to_bits() != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
        { return Err(invalid("aperture performance requires a matching clock, positive capacity/window and an unadvanced complete physical system")); }
        let exterior_gas = match (&system, observation) {
            (CoupledAperture::Dg(model), ApertureObservation::TubeBaffled(_)) => {
                let record = model.load().radiation_record()
                    .ok_or_else(|| invalid("DG exterior observation requires its computed baffled reaction source"))?;
                if record.apertures().is_empty()
                    || record.apertures().len() != model.load().terminal_ports().len()
                {
                    return Err(invalid("DG exterior observation requires every original bound terminal opening"));
                }
                let bound = *record.gas_record().state();
                if exterior_gas.is_some_and(|gas| gas != bound) {
                    return Err(invalid("DG receiver exterior differs from its computed reaction source"));
                }
                Some(bound)
            }
            _ => exterior_gas,
        };
        if let Some(exterior) = exterior_gas {
            exterior.validate_coherence().map_err(|e|
                GestureCompileError::Render(receiver_error(e.to_string())))?;
            match (&system, observation) {
                (CoupledAperture::Network(network), ApertureObservation::NetworkBaffled {node, ..}) => {
                    let gases = network.section_gases()
                        .ok_or_else(|| invalid("independent exterior requires complete internal gas states"))?;
                    let mut adjacent = network.spec().sections.iter().enumerate()
                        .filter(|(_, s)| s.nodes.contains(&node));
                    let (index, _) = adjacent.next()
                        .ok_or_else(|| invalid("exterior observation has no adjacent section"))?;
                    if adjacent.next().is_some() || gases[index].pressure != exterior.pressure {
                        return Err(invalid("exterior observation requires a degree-one port and common static pressure; mean flow is not modeled"));
                    }
                }
                (CoupledAperture::Dg(model), ApertureObservation::TubeBaffled(_)) => {
                    if model.load().duct().outlet_gas().pressure != exterior.pressure {
                        return Err(invalid("DG receiver and outlet must share static pressure"));
                    }
                }
                _ => return Err(invalid("independent exterior requires a bound outlet observation")),
            }
        }
        let receiver = match (&system, observation) {
            (_, ApertureObservation::Inlet) | (CoupledAperture::Tube(_), ApertureObservation::TubeTerminal) => None,
            (CoupledAperture::Network(m), ApertureObservation::NetworkNode(n)) if n < m.spec().nodes.len() => None,
            (CoupledAperture::Dg(model), ApertureObservation::TubeBaffled(location)) => {
                Some(model.load().baffled_receiver(config.sample_rate_hz, location)
                    .map_err(|e| GestureCompileError::Render(RenderError::Voice(e)))?)
            }
            (CoupledAperture::Tube(m), ApertureObservation::TubeBaffled(location)) => {
                // Radius and gas are the actual physical tube's, never an
                // independently tuned observer area or invented gain. Waveguide
                // construction is quiescent, so the initial flow history is zero.
                Some(BaffledPressure::circular_outlet(m.spec().radius_m,
                    config.sample_rate_hz, location, RayleighMedium {
                        density: a.spec().density_kg_m3, sound_speed: m.spec().sound_speed_m_s,
                    }).map_err(|e| GestureCompileError::Render(receiver_error(e)))?)
            }
            (CoupledAperture::Network(m), ApertureObservation::NetworkBaffled {node,receiver}) => {
                use super::network::NetworkNode;
                if !matches!(m.spec().nodes.get(node), Some(NetworkNode::Termination {..}
                    | NetworkNode::Impedance {..} | NetworkNode::Relaxation {..})) {
                    return Err(invalid("baffled network observation requires a degree-one passive terminal"));
                }
                let (index, section)=m.spec().sections.iter().enumerate().find(|(_,s)|s.nodes.contains(&node))
                    .ok_or_else(||invalid("baffled network terminal has no physical section"))?;
                let (density, sound_speed)=m.section_medium(index).expect("existing outlet section");
                let (density, sound_speed) = exterior_gas.map_or((density, sound_speed),
                    |gas| (gas.density, gas.sound_speed));
                // Omission is the exact original local-gas choice. An explicit
                // exterior describes the pressure-matched outlet environment,
                // not a replacement for the bore's propagation properties.
                Some(BaffledPressure::circular_outlet(section.radius_m,config.sample_rate_hz,receiver,
                    RayleighMedium {density,sound_speed})
                    .map_err(|e|GestureCompileError::Render(receiver_error(e)))?)
            }
            _ => return Err(invalid("pressure observation is not a node or terminal of the supplied system")),
        };
        let controls = compile_schedule(&schedule, &config)?;
        let dg_terminal_flows = match (&system, &receiver) {
            (CoupledAperture::Dg(model), Some(_)) => vec![0.0; model.load().terminal_ports().len()],
            _ => Vec::new(),
        };
        Ok(Self { system, observation, receiver, dg_terminal_flows, exterior_gas, config, schedule, controls, next: 0,
            held_pressure_pa: 0.0, poisoned: false, forces: None, boundaries: None,
            boundary_frame: None, mechanical_work_j: 0.0,
            finite_gas: None,
            finite_transport: None,
            finite_wall: None,
            observation_policy: ObservationPolicy::Diagnostic })
    }

    /// Read-only physical system; callers cannot step it behind the sample clock.
    #[must_use]
    pub const fn system(&self) -> &CoupledAperture { &self.system }
    /// Original source gestures, unchanged by audio resampling or block partition.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Exact finite clock/compilation choices.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Named physical location of the output trace.
    #[must_use]
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Immutable receiver geometry/delay information, absent for internal traces.
    /// The physical propagation delay remains inside this source; downstream
    /// rate-conversion or ensemble alignment must not remove it.
    #[must_use]
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Explicitly bound exterior, absent on the legacy local-gas path.
    #[must_use]
    pub const fn exterior_gas(&self) -> Option<&GasState> { self.exterior_gas.as_ref() }
    /// Actual density and sound speed used by the receiver, not the source inlet.
    #[must_use]
    pub fn baffled_medium(&self) -> Option<RayleighMedium> {
        if self.receiver.is_none() { return None; }
        if let Some(gas) = self.exterior_gas {
            return Some(RayleighMedium { density: gas.density, sound_speed: gas.sound_speed });
        }
        match (&self.system, self.observation) {
            (CoupledAperture::Tube(t), ApertureObservation::TubeBaffled(_)) =>
                Some(RayleighMedium { density: t.aperture().spec().density_kg_m3,
                    sound_speed: t.spec().sound_speed_m_s }),
            (CoupledAperture::Network(n), ApertureObservation::NetworkBaffled {node, ..}) => {
                let (index, _) = n.spec().sections.iter().enumerate()
                    .find(|(_, s)| s.nodes.contains(&node))?;
                let (density, sound_speed) = n.section_medium(index)?;
                Some(RayleighMedium { density, sound_speed })
            }
            _ => None,
        }
    }
    /// The complete gas behind [`Self::baffled_medium`], when the receiver hears a
    /// sourced state: the bound exterior, or the observed network section's own
    /// gas. The legacy tube spec and a coefficient-only network carry none.
    #[must_use]
    pub fn baffled_gas(&self) -> Option<&GasState> {
        self.receiver.as_ref()?;
        if let Some(gas) = self.exterior_gas.as_ref() { return Some(gas); }
        match (&self.system, self.observation) {
            (CoupledAperture::Network(n), ApertureObservation::NetworkBaffled {node, ..}) => {
                let index = n.spec().sections.iter().position(|s| s.nodes.contains(&node))?;
                n.section_gases()?.get(index)
            }
            _ => None,
        }
    }
    /// Controls committed with successful mechanical samples, never ahead of them.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Remaining controls on the original mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Finite remaining source samples, not resampled output frames.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.samples_rendered() }

    fn step_pressure(&mut self, pressure_pa: f64, force_n: f64)
        -> Result<(f64, f64, Option<BoundaryNetworkFrame>), RenderError> {
        let drive = TubeDrive { upstream_pressure_pa: pressure_pa, body_flow_m3_s: 0.0 };
        let boundary_drives = self.boundaries.as_ref().map_or(&[][..], |b| b.candidate());
        match &mut self.system {
            CoupledAperture::Dg(m) => {
                let frame = if let Some(transport) = &mut self.finite_transport {
                    transport.step(m, drive, force_n)?
                } else if let Some(gas) = &mut self.finite_gas {
                    gas.step(m, drive, force_n)?
                } else if let Some(wall) = &mut self.finite_wall {
                    wall.step(m, drive, force_n)?
                } else {
                    m.step_with_force(drive, force_n).map_err(RenderError::Voice)?
                };
                let pressure = match self.observation {
                    ApertureObservation::Inlet => frame.aperture.bore_pressure_pa,
                    ApertureObservation::TubeBaffled(_) => {
                        if frame.terminals.len() != self.dg_terminal_flows.len()
                            || frame.terminals.is_empty()
                        {
                            return Err(sizing("DG microphone requires every accepted terminal flow"));
                        }
                        for (flow, terminal) in self.dg_terminal_flows.iter_mut().zip(&frame.terminals) {
                            *flow = terminal.outward_flow_m3_s;
                        }
                        self.receiver.as_mut().expect("admitted DG outlet receiver")
                            .step(&self.dg_terminal_flows).map_err(receiver_error)?
                    }
                    _ => unreachable!("observer admitted against the exclusively owned system"),
                };
                Ok((pressure, frame.aperture.mechanical_work_j, None))
            }
            CoupledAperture::Tube(m) => {
                let f = m.step_with_force(drive, force_n).map_err(RenderError::Voice)?;
                let pressure = match self.observation {
                    ApertureObservation::Inlet => f.aperture.bore_pressure_pa,
                    ApertureObservation::TubeTerminal => f.waveguide.terminal_pressure_pa,
                    ApertureObservation::TubeBaffled(_) => self.receiver.as_mut()
                        .expect("baffled observation constructed its complete receiver")
                        .step(&[f.waveguide.terminal_flow_m3_s]).map_err(receiver_error)?,
                    _ => unreachable!("observer admitted against the exclusively owned system"),
                };
                Ok((pressure, f.aperture.mechanical_work_j, None))
            }
            CoupledAperture::Network(m) => {
                let f = m.step_with_boundaries(drive, force_n, boundary_drives).map_err(RenderError::Voice)?;
                let pressure = match self.observation {
                    ApertureObservation::Inlet => f.aperture.bore_pressure_pa,
                    ApertureObservation::NetworkNode(n) => m.node_frame(n)
                        .expect("admitted node remains in the owned topology").pressure_pa,
                    ApertureObservation::NetworkBaffled {node,..} => {
                        let flow=m.node_frame(node).expect("admitted terminal remains in owned topology")
                            .net_flow_into_node_m3_s;
                        self.receiver.as_mut().expect("admitted network outlet receiver")
                            .step(&[flow]).map_err(receiver_error)?
                    }
                    _ => unreachable!("observer admitted against the exclusively owned system"),
                };
                Ok((pressure, f.aperture.mechanical_work_j, f.boundaries))
            }
        }
    }
}
impl PressureRenderer for AperturePerformance {
    fn samples_rendered(&self) -> u64 { self.system.aperture().accepted_steps() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz { return Err(sizing("pressure output must use the physical aperture clock")); }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if matches!(&self.system, CoupledAperture::Network(m) if m.boundary_apertures().is_some())
            && self.boundaries.is_none() {
            return Err(sizing("keyed aperture performance requires explicit boundary controls before rendering"));
        }
        if samples > self.remaining_samples() { return Err(sizing("pressure request exceeds the complete aperture performance window")); }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("aperture callback exceeds its construction-time capacity")); }
        for slot in output {
            let sample = self.samples_rendered();
            let (next, pressure) = pressure_candidate(&self.controls, self.next, self.held_pressure_pa, sample);
            let force_n = self.forces.as_ref().map_or(0.0, |f| f.candidate(sample));
            if let Some(boundaries) = &mut self.boundaries { boundaries.prepare(sample); }
            match self.step_pressure(pressure, force_n) {
                Ok((value, work, boundary_frame)) => {
                    *slot = value; self.next = next; self.held_pressure_pa = pressure;
                    if let Some(forces) = &mut self.forces { forces.accept(sample); }
                    if let Some(boundaries) = &mut self.boundaries { boundaries.accept(); }
                    if let Some(gas) = &mut self.finite_gas { gas.accept(sample); }
                    if let Some(transport) = &mut self.finite_transport { transport.accept(sample); }
                    if let Some(wall) = &mut self.finite_wall { wall.accept(sample); }
                    self.mechanical_work_j = work;
                    self.boundary_frame = boundary_frame;
                }
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}

// Shared source-clock admission/compilation for scalar and independent blades.
// Physical owners remain different; control semantics and work budgets do not.
pub(in crate::bernoulli_aperture) fn compile_schedule(schedule: &GestureSchedule, config: &AperturePerformanceConfig)
    -> Result<Vec<ScheduledControl>, GestureCompileError>
{
    let [track] = schedule.tracks() else {
        return Err(invalid("one explicitly typed blowing-pressure track is required"));
    };
    let rate = schedule.control_rate_hz;
    if rate == 0 || rate > config.sample_rate_hz {
        return Err(invalid("control rate must be positive and no faster than mechanics"));
    }
    let last_tick = u128::from(config.samples - 1) * u128::from(rate)
        / u128::from(config.sample_rate_hz);
    if last_tick + 1 > config.max_controls as u128 {
        return Err(invalid("pressure sampling exceeds the declared control-tick/storage budget"));
    }
    if track.events.iter().any(|e| e.time_s > last_tick as f64 / f64::from(rate)) {
        return Err(invalid("pressure command begins after the last observed control tick"));
    }
    let controls=compile_pressure_gestures(schedule, &[PressureGestureBinding {
        track: track.id.clone(), voice: 0,
    }], config.sample_rate_hz, config.samples, config.max_compile_work)?;
    if controls.len()>config.max_controls {
        return Err(invalid("compiled pressure exceeds the declared control storage allowance"));
    }
    Ok(controls)
}

// Preview controls without publishing them ahead of their physical sample.
pub(in crate::bernoulli_aperture) fn pressure_candidate(controls: &[ScheduledControl], mut next: usize, mut pressure: f64,
    sample: u64) -> (usize, f64)
{
    while let Some(event) = controls.get(next) {
        if event.sample != sample { break; }
        if let ControlDelta::SetBlowingPressure { pressure_pa, .. } = event.delta {
            pressure = pressure_pa;
        } else { unreachable!("existing compiler emits only pressure assignments"); }
        next += 1;
    }
    (next, pressure)
}

#[cfg(test)]
mod exterior_tests;

#[cfg(test)]
mod qualification;

#[cfg(test)]
mod dg_tests;
