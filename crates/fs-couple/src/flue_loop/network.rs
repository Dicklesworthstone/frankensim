//! A physical tube graph under the existing jet-drive island.
//!
//! No second jet equation or propagation integrator. Original node addresses
//! survive viscothermal subdivision, so loads, branches and outlet observations
//! refer to the same geometry. Gas fields are frozen, common-pressure states;
//! mixing, mean-flow convection and time-varying materials are not modeled.

use super::{FlueError, FlueGeometry, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{SectionRealization, TubeNetworkSpec};
use crate::bernoulli_aperture::network::boundaries::{
    BoundaryAperture, BoundaryApertures, BoundaryDrive, BoundaryNetworkFrame,
};
use crate::bernoulli_aperture::fractional::{FractionalSection, FractionalTubePolicy};
use crate::bernoulli_aperture::viscothermal::{
    ViscothermalSection, ViscothermalSectionSpec, with_regional_viscothermal_sections,
    with_regional_fractional_viscothermal_sections,
};
use crate::bernoulli_aperture::wall::WallState;
use crate::bernoulli_aperture::wall::material::{
    BoundMaterialWall, MaterialWallInput, bind_material_walls,
};
use crate::bernoulli_aperture::radiation::environment::multiport::{
    MimoEnvironment, MimoRadiationLimits, MimoRadiationModel,
};
use crate::wind_line::WindLineBank;
use fs_aeroac::jetcard::JetCard;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_scenario::assembly::{ParameterKey, ParameterPath};
use fs_scenario::entity::EntityId;
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
    /// Declared minimum jet speed U_min [m/s]; it sizes the retained convection history. Below
    /// it, a jet whose lookup would read motion older than that history carries no deflection,
    /// and the windway's steady flow splits at zero deflection (ruling 8132). A validity edge,
    /// not a stability limit. Zero blowing keeps its own path and preserves stored vibration.
    pub minimum_jet_speed_m_s: f64,
    /// Explicit incompressible-jet Mach ceiling, in (0, 0.3].
    pub maximum_jet_mach: f64,
    /// Maximum retained jet convection samples, including interpolation slots.
    /// Separate from graph payload; bounded above by 1,048,576 samples.
    pub max_jet_history_samples: usize,
}

/// An explicitly actuated physical boundary on the same source gas and graph.
/// `closing_force` is a dimensioned description path, never a fingering preset.
pub struct FlueBoundaryAperture {
    /// Original driven endpoint and existing moving-slit mechanics.
    pub boundary: BoundaryAperture,
    /// Entity/parameter address of the nonnegative closing force [N].
    pub closing_force: ParameterPath,
    /// Initial externally applied closing force [N].
    pub initial_closing_force_n: f64,
    /// External reservoir pressure relative to the common equilibrium [Pa].
    pub exterior_pressure_pa: f64,
}

/// Exact typed source controls consumed by the existing pressure/angle owner.
/// A v2 performance contains only resolved paths; it never mixes legacy targets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlueSourcePaths {
    /// Source blowing pressure in pascals.
    pub pressure: ParameterPath,
    /// Source jet direction in radians.
    pub jet_angle: ParameterPath,
}
impl FlueSourcePaths {
    /// Bind both physical parameters to an explicitly declared source entity.
    #[must_use]
    pub fn for_entity(entity: EntityId) -> Self {
        Self {
            pressure: ParameterPath { entity, parameter: ParameterKey::new("blowing_pressure", fs_scenario::gesture::path::PRESSURE) },
            jet_angle: ParameterPath { entity, parameter: ParameterKey::new("jet_angle", fs_qty::Dims::NONE) },
        }
    }
}

/// The existing wave/load owner and its exact original material/geometry inputs.
/// The jet card's estimate-only authority is not changed by this construction.
pub struct FlueNetwork {
    source: FlueNetworkSpec,
    wall_input: Option<MaterialWallInput>,
    walls: Vec<BoundMaterialWall>,
    graph: TubeNetworkSpec,
    gases: Vec<GasState>,
    represented: Vec<SectionRealization>,
    losses: Vec<ViscothermalSection>,
    physics: WaveguideNetwork,
    last: Option<NetworkFrame>,
    fractional: Option<(FractionalTubePolicy, Vec<FractionalSection>)>,
    mimo_environment: Option<MimoEnvironment>,
    boundaries: Option<BoundaryApertures>,
    boundary_paths: Vec<ParameterPath>,
    boundary_drives: Vec<BoundaryDrive>,
    source_paths: Option<FlueSourcePaths>,
    opening_losses: Vec<crate::bernoulli_aperture::opening_loss::BoundOpeningLoss>,
}

pub(super) enum FlueResonator {
    Tmm(WindLineBank),
    Network(Box<FlueNetwork>),
    Dg(Box<super::dg::FlueDg>),
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
    fn new(source: FlueNetworkSpec, wall_input: Option<MaterialWallInput>, dt: f64,
        zc: f64, gate: &CancelGate, policy: Option<FractionalTubePolicy>) -> Result<Self, FlueError>
    {
        poll(gate)?;
        // Bind to ORIGINAL geometry before subdivision. Both operations preserve
        // original node indices; material area must not depend on loss-cell count.
        let (wall_graph, walls) = match &wall_input {
            Some(input) => {
                // Bound duplicate source storage BEFORE cloning the receipts.
                if input.max_walls > 1024 || input.assignments.len() > input.max_walls {
                    return Err(invalid("material wall source storage exceeds its assignment budget"));
                }
                if input.band_hz[1] >= 0.5 / dt {
                    return Err(invalid("material wall applicability band must lie below the source Nyquist frequency"));
                }
                bind_material_walls(source.graph.clone(), input.assignments.clone(),
                    input.band_hz, input.max_walls, gate).map_err(physical)?
            }
            None => (source.graph.clone(), Vec::new()),
        };
        for wall in &walls { poll(gate)?; wall.validate_time_step(dt).map_err(physical)?; }
        let (graph, gases, losses) = match policy {
            Some(policy) => with_regional_fractional_viscothermal_sections(
                wall_graph, source.gases.clone(), dt, &source.viscothermal, policy, gate),
            None => with_regional_viscothermal_sections(
                wall_graph, source.gases.clone(), dt, &source.viscothermal, gate),
        }.map_err(physical)?;
        // A caller cannot provide a second, inconsistent source cross-section.
        let actual = graph.inlet_impedance_with_gases(&gases).map_err(physical)?;
        if actual.to_bits() != zc.to_bits() {
            return Err(invalid("flue mouth bore area must match the physical inlet section"));
        }
        let (physics, represented, fractional) = if let Some(policy) = policy {
            let (physics, represented, details) = graph
                .realize_fractional_waveguide(&gases, dt, policy, gate).map_err(physical)?;
            (physics, represented, Some((policy, details)))
        } else {
            let (physics, represented) = graph.realize_waveguide(&gases, dt, gate).map_err(physical)?;
            (physics, represented, None)
        };
        poll(gate)?;
        Ok(Self { source, wall_input, walls, graph, gases, represented, losses, physics,
            last: None, fractional, mimo_environment: None, boundaries: None,
            boundary_paths: Vec::new(), boundary_drives: Vec::new(), source_paths: None,
            opening_losses: Vec::new() })
    }
    /// Original data before any loss-cell subdivision.
    #[must_use]
    pub const fn source(&self) -> &FlueNetworkSpec { &self.source }
    /// Original solid assignments, frequency applicability and work limit.
    /// Absent for legacy/raw-coefficient construction; never invented from loads.
    #[must_use]
    pub const fn wall_input(&self) -> Option<&MaterialWallInput> { self.wall_input.as_ref() }
    /// Original material states and their geometry-derived pressure/flow laws.
    #[must_use]
    pub fn material_walls(&self) -> &[BoundMaterialWall] { &self.walls }
    /// Nonlinear terminal losses with their original local gas and opening area.
    pub fn opening_losses(&self) -> &[crate::bernoulli_aperture::opening_loss::BoundOpeningLoss] {
        &self.opening_losses
    }
    /// Accepted displacement, velocity and BASE mechanical storage at a material node.
    /// Hereditary spring energy is separately available via `wall_relaxation_energy_j`.
    /// Reads the same shunt advanced by the jet/bore solve, not another oscillator.
    ///
    /// # Errors
    /// Unbound node, mismatching node law, or nonfinite motion.
    pub fn wall_state(&self, node: usize) -> Result<WallState, FlueError> {
        let wall = self.walls.iter().find(|wall| wall.source().node == node)
            .ok_or_else(|| invalid("wall motion requires a source-bound material node"))?;
        wall.observe_waveguide(&self.physics).map_err(physical)
    }
    /// Accepted hereditary spring storage [J], already included in network energy.
    /// Adding this to `wall_state(node)?.stored_energy_j` gives complete wall storage.
    ///
    /// # Errors
    /// Unbound node, mismatched material load or unavailable/nonfinite memory.
    pub fn wall_relaxation_energy_j(&self, node: usize) -> Result<f64, FlueError> {
        let wall = self.walls.iter().find(|wall| wall.source().node == node)
            .ok_or_else(|| invalid("wall memory requires a source-bound material node"))?;
        wall.relaxation_stored_energy_j(&self.physics).map_err(physical)
    }
    /// Actual lowered graph. Original node indices remain valid.
    #[must_use]
    pub const fn graph(&self) -> &TubeNetworkSpec { &self.graph }
    /// Complete ordered states actually used for propagation and observation.
    #[must_use]
    pub fn section_gases(&self) -> &[GasState] { &self.gases }
    /// Physical lengths and impedances admitted by the shared geometry owner.
    #[must_use]
    pub fn represented_sections(&self) -> &[SectionRealization] { &self.represented }
    /// Selected phase policy and residual transit of each actually lowered section.
    /// None is the original integer image. Selected viscothermal cells preserve
    /// original physical lengths; their separate loss/error records remain visible.
    #[must_use]
    pub fn fractional_sections(&self) -> Option<(FractionalTubePolicy, &[FractionalSection])> {
        self.fractional.as_ref().map(|(policy, details)|(*policy, details.as_slice()))
    }
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
    pub fn stored_energy_j(&self) -> f64 {
        self.physics.stored_energy_j() + self.boundaries.as_ref().map_or(0.0, BoundaryApertures::stored_energy_j)
    }
    /// Original boundary mechanics, contact and retained material histories.
    pub const fn boundary_apertures(&self) -> Option<&BoundaryApertures> { self.boundaries.as_ref() }
    /// Explicitly bound pressure/angle consumers for all-typed v2 performances.
    pub const fn source_control_paths(&self) -> Option<&FlueSourcePaths> { self.source_paths.as_ref() }
    /// Physical node addressed by an explicitly bound closing-force parameter.
    pub fn boundary_node(&self, path: &ParameterPath) -> Option<usize> {
        self.boundary_paths.iter().position(|p| p == path).map(|i| self.boundary_drives[i].node)
    }
    /// Current prescribed reservoir pressure and signed mechanical force. A
    /// control changes this input only; state advances with the shared next step.
    pub fn boundary_drive(&self, node: usize) -> Option<&BoundaryDrive> {
        self.boundary_drives.iter().find(|drive| drive.node == node)
    }
    /// Last combined network/boundary work, including independent exterior loads.
    pub fn boundary_frame(&self) -> Option<&BoundaryNetworkFrame> {
        self.boundaries.as_ref()?.last_frame()
    }
    /// Source-bound mutual terminal operator actually installed in this network.
    #[must_use]
    pub const fn mimo_environment(&self) -> Option<&MimoEnvironment> { self.mimo_environment.as_ref() }
    /// Same-sample pressure, flow and signed work of the shared exterior model.
    #[must_use]
    pub fn mimo_port_frames(&self) -> Option<&[fs_vfit::waveguide::network::multiport::MimoPortFrame]> {
        self.physics.mimo_port_frames()
    }
    /// Continuous certificate and separate actual-ZOH sampled-response diagnostic.
    #[must_use]
    pub fn mimo_passivity(&self) -> Option<(&fs_vfit::mimo::MimoPassivityReport, &fs_vfit::mimo::MimoPassivityReport)> {
        self.physics.mimo_passivity()
    }
}

impl FlueResonator {
    pub(super) fn incoming(&self) -> f64 {
        match self {
            Self::Tmm(line) => line.incoming(),
            Self::Network(n) => n.physics.incoming_pressure_pa(),
            Self::Dg(duct) => duct.incoming(),
        }
    }
    pub(super) fn push(&mut self, mouth: super::MouthStep) -> Result<f64, FlueError> {
        let outgoing = mouth.outgoing_pa;
        match self {
            Self::Tmm(line) => line.push(outgoing).map_err(FlueError::Line),
            Self::Dg(duct) => duct.push(mouth),
            Self::Network(n) => {
                let frame = if let Some(boundaries) = &mut n.boundaries {
                    boundaries.step(&mut n.physics, outgoing, &n.boundary_drives)
                        .map_err(physical)?.network
                } else { n.physics.step(outgoing).map_err(wave)? };
                n.last = Some(frame);
                Ok(n.physics.incoming_pressure_pa())
            }
        }
    }
    pub(super) fn switch_fingering(&mut self, index: usize, fade: usize) -> Result<(), FlueError> {
        match self {
            Self::Tmm(line) => line.switch_fingering(index, fade).map_err(FlueError::Line),
            Self::Network(_) | Self::Dg(_) => Err(invalid("a physical flue bore has no TMM fingering bank")),
        }
    }
    // A positive jet below the declared minimum speed is admitted: the voice reads the rest state
    // or takes the declared state there (ruling 8132).
    pub(super) fn validate_jet_speed(&self, speed: f64) -> Result<(), FlueError> {
        if let Self::Dg(duct) = self { return duct.validate_jet_speed(speed); }
        if let Self::Network(n) = self {
            let source = &n.source;
            let inlet = source.graph.inlet_section_index().map_err(physical)?;
            if !speed.is_finite()
                || speed > source.maximum_jet_mach * source.gases[inlet].sound_speed {
                return Err(invalid("blowing pressure exceeds the declared low-Mach jet domain"));
            }
        }
        Ok(())
    }
}

impl FlueVoice {
    /// Attach the actual source's entity paths before rendering. The legacy v1
    /// targets remain available; every v2 track must resolve an explicit path.
    pub fn with_source_control_paths(mut self, paths: FlueSourcePaths) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0
            || paths != FlueSourcePaths::for_entity(paths.pressure.entity) {
            return Err(invalid("source controls require cold, coherent pressure/angle parameter paths"));
        }
        let FlueResonator::Network(network) = &mut self.line else {
            return Err(invalid("typed source controls require the explicit physical graph"));
        };
        if network.source_paths.is_some() { return Err(invalid("source control paths are already bound")); }
        network.source_paths = Some(paths);
        Ok(self)
    }

    /// Attach material/contact boundaries before the first source step. Each
    /// opening keeps its full original chimney and mechanical histories during
    /// later force changes. This exposes the existing moving-slit approximation;
    /// it does not infer circular pads, key levers or their material coefficients.
    pub fn with_boundary_apertures(mut self, inputs: Vec<FlueBoundaryAperture>,
        max_apertures: usize, gate: &CancelGate) -> Result<Self, FlueError>
    {
        poll(gate)?;
        if self.poisoned || self.samples != 0 || !(1..=64).contains(&max_apertures)
            || inputs.is_empty() || inputs.len() > max_apertures {
            return Err(invalid("moving boundaries require a cold flue and an explicit count allowance"));
        }
        let FlueResonator::Network(network) = &mut self.line else {
            return Err(invalid("moving boundaries require the explicit physical tube graph"));
        };
        if network.boundaries.is_some() { return Err(invalid("moving boundaries are already bound")); }
        let mut paths = Vec::with_capacity(inputs.len());
        let mut drives = Vec::with_capacity(inputs.len());
        let mut mechanics = Vec::with_capacity(inputs.len());
        for input in inputs {
            poll(gate)?;
            if input.closing_force.parameter.name != "closing_force"
                || input.closing_force.parameter.dims != fs_scenario::gesture::path::FORCE
                || paths.contains(&input.closing_force)
                || !input.initial_closing_force_n.is_finite() || input.initial_closing_force_n < 0.0
                || !input.exterior_pressure_pa.is_finite() {
                return Err(invalid("moving boundaries need distinct force-dimensioned closing_force paths and finite physical drives"));
            }
            paths.push(input.closing_force);
            drives.push(BoundaryDrive { node: input.boundary.node,
                exterior_pressure_pa: input.exterior_pressure_pa,
                opening_force_n: -input.initial_closing_force_n });
            mechanics.push(input.boundary);
        }
        let boundaries = BoundaryApertures::new(&network.graph, &network.gases,
            &network.physics, mechanics, max_apertures).map_err(physical)?;
        poll(gate)?;
        network.boundaries = Some(boundaries);
        network.boundary_paths = paths;
        network.boundary_drives = drives;
        Ok(self)
    }

    /// Admission of a source closing-force assignment without changing state.
    pub fn validate_boundary_closing_force(&self, node: usize, force_n: f64) -> Result<(), FlueError> {
        let FlueResonator::Network(network) = &self.line else {
            return Err(invalid("boundary force requires an explicitly bound physical graph"));
        };
        if !force_n.is_finite() || force_n < 0.0
            || !network.boundary_drives.iter().any(|drive| drive.node == node) {
            return Err(invalid("boundary closing force must be nonnegative finite newtons at a bound node"));
        }
        Ok(())
    }

    pub(super) fn set_boundary_closing_force(&mut self, node: usize, force_n: f64) -> Result<(), FlueError> {
        self.validate_boundary_closing_force(node, force_n)?;
        let FlueResonator::Network(network) = &mut self.line else { unreachable!("admitted graph"); };
        network.boundary_drives.iter_mut().find(|drive| drive.node == node)
            .expect("admitted boundary").opening_force_n = -force_n;
        Ok(())
    }

    /// Attach a source-bound mutual radiation operator to an unadvanced physical
    /// jet/bore, preserving its original walls, gases, losses, windway and seed.
    /// One matrix solve returns pressure to every selected outlet in the SAME
    /// sample; all source jets continue through the existing flue runtime.
    ///
    /// The declaration's exterior and aperture order are retained for every
    /// later FlueOutletPerformance receiver. Changing a baked material requires
    /// its corresponding source operator, never an arbitrary gain rescaling.
    /// The actual ZOH grid is explicitly GridOnly: no continuous certificate is
    /// promoted to a discrete storage law or complete M2.5-W3c qualification.
    /// Full clock bandwidth is mandatory. Delayed/improper operators, existing
    /// reactive loads, live attachment and non-baffled observer geometry refuse.
    pub fn with_mimo_environment(mut self, source: MimoRadiationModel,
        radiation: crate::pcm_wav::baffled::openings::duct::DuctRadiationSpec,
        limits: MimoRadiationLimits, gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        poll(gate)?;
        if self.poisoned || self.samples != 0
            || self.dt.to_bits() != (1.0 / f64::from(radiation.reconstruction.mechanics_rate_hz())).to_bits() {
            return Err(invalid("MIMO exterior requires an unadvanced flue and its actual mechanical clock"));
        }
        let FlueResonator::Network(network) = &mut self.line else {
            return Err(invalid("mutual terminal radiation requires an explicitly ported physical tube graph"));
        };
        let (environment, terminal) = MimoEnvironment::prepare(&network.graph, &network.gases,
            source, radiation, limits, gate).map_err(physical)?;
        network.physics.install_mimo_terminals(terminal).map_err(wave)?;
        poll(gate)?;
        network.mimo_environment = Some(environment);
        Ok(self)
    }

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
        Self::from_network_with_walls(card, displacement_gain, seed_amplitude_rel, geometry,
            source, None, sample_rate_hz, voice_seed, voice_index, gate, None)
    }

    /// The same jet and regional gas network with source-bound structural walls.
    ///
    /// Wall displacement returns volume flow into the SAME pressure junction that
    /// drives the jet. Original geometry/material receipts survive viscothermal
    /// subdivision, source stepping and `FlueOutletPerformance` rendering. Changing
    /// a wall card changes mechanical mass/stiffness/loss, never a fitted sound gain.
    /// The no-wall constructor retains its original path and public signature.
    ///
    /// Only the thin, small-motion isotropic membrane image in `wall::material`
    /// is admitted. No axial shell propagation, orthotropic wood, dynamic material
    /// change, exterior shell radiation, or measured flute-fidelity claim follows.
    /// The material band is not an antialiasing filter; nonlinear out-of-band
    /// content and displacement validity still need application-level assessment.
    ///
    /// # Errors
    /// Original jet/network refusals, missing or incompatible wall data/geometry,
    /// a wall band at/above Nyquist, exhausted admission limits or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_network(
        card: Option<&JetCard>, displacement_gain: f64, seed_amplitude_rel: f64,
        geometry: FlueGeometry, source: FlueNetworkSpec, walls: MaterialWallInput,
        sample_rate_hz: u32, voice_seed: u64, voice_index: u32, gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        Self::from_network_with_walls(card, displacement_gain, seed_amplitude_rel, geometry,
            source, Some(walls), sample_rate_hz, voice_seed, voice_index, gate, None)
    }

    /// Keep sub-sample bore travel in the ORIGINAL jet-drive and material-wall
    /// runtime. The supplied gas determines both source flow and each section's
    /// L/c delay; no retuned oscillator, additional jet law or receiver gain.
    /// The existing outlet performers consume this same voice and accepted flows.
    ///
    /// Both wall assignment and transit selection happen before any sample. The
    /// original integer constructors remain unchanged. A fixed first-order
    /// allpass has finite-band phase error and stored memory, not exact broadband
    /// delay. It cannot be changed live without a separate boundary-work law.
    ///
    /// Selected viscothermal cells use original physical lengths and regional
    /// gas properties through the shared fractional loss lowering. The phase
    /// band must cover each selected loss band; both the phase budget and the
    /// original physical/scattering gates apply. No selected loss or structural
    /// wall is silently removed, and material wall area precedes subdivision.
    ///
    /// # Errors
    /// Original jet/material/geometry refusals or fractional phase/memory limits.
    #[allow(clippy::too_many_arguments)]
    pub fn from_fractional_network(
        card: Option<&JetCard>, displacement_gain: f64, seed_amplitude_rel: f64,
        geometry: FlueGeometry, source: FlueNetworkSpec, walls: Option<MaterialWallInput>,
        policy: FractionalTubePolicy, sample_rate_hz: u32,
        voice_seed: u64, voice_index: u32, gate: &CancelGate,
    ) -> Result<Self, FlueError> {
        Self::from_network_with_walls(card, displacement_gain, seed_amplitude_rel, geometry,
            source, walls, sample_rate_hz, voice_seed, voice_index, gate, Some(policy))
    }

    #[allow(clippy::too_many_arguments)]
    fn from_network_with_walls(
        card: Option<&JetCard>, displacement_gain: f64, seed_amplitude_rel: f64,
        geometry: FlueGeometry, source: FlueNetworkSpec, walls: Option<MaterialWallInput>,
        sample_rate_hz: u32, voice_seed: u64, voice_index: u32, gate: &CancelGate,
        policy: Option<FractionalTubePolicy>,
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
            |zc| FlueNetwork::new(source, walls, 1.0 / f64::from(sample_rate_hz), zc, gate, policy)
                .map(|n| FlueResonator::Network(Box::new(n))))?;
        poll(gate)?;
        Ok(voice)
    }
    /// Read the actual regional tube state; absent for the legacy TMM image.
    #[must_use]
    pub fn network(&self) -> Option<&FlueNetwork> {
        match &self.line {
            FlueResonator::Network(n) => Some(n),
            FlueResonator::Tmm(_) | FlueResonator::Dg(_) => None,
        }
    }
    /// Add explicit gas/geometry-derived separation losses at physical outlets.
    /// Original radiation, wall and wave histories share the same solved flow.
    /// This is a cold quasi-steady closure, not resolved vortex shedding.
    pub fn with_opening_losses(mut self,
        inputs: Vec<crate::bernoulli_aperture::opening_loss::OpeningLossInput>,
        gate: &CancelGate) -> Result<Self, FlueError> {
        if self.samples != 0 { return Err(invalid("opening losses require an unadvanced flue source")); }
        let FlueResonator::Network(network) = &mut self.line else {
            return Err(invalid("opening losses require the retained physical tube graph"));
        };
        if !network.opening_losses.is_empty() {
            return Err(invalid("opening losses must be supplied together once"));
        }
        network.opening_losses = crate::bernoulli_aperture::opening_loss::bind_opening_losses(
            &network.graph, &network.gases, &mut network.physics, inputs, gate).map_err(physical)?;
        Ok(self)
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
            FlueResonator::Tmm(_) | FlueResonator::Dg(_) => Err(invalid("terminal controls require the physical graph image")),
        }
    }
    /// Validate a held physical pressure before compiling a source schedule.
    pub fn validate_blowing_pressure(&self, pressure_pa: f64) -> Result<(), FlueError> {
        self.validated_jet_speed(pressure_pa).map(|_| ())
    }
    pub(super) fn validated_jet_speed(&self, pressure_pa: f64) -> Result<f64, FlueError> {
        if self.poisoned { return Err(invalid("flue voice is poisoned by a previous failed sample")); }
        if !pressure_pa.is_finite() || pressure_pa < 0.0 {
            return Err(invalid("blowing pressure must be finite and non-negative"));
        }
        let speed = self.speed_at_pressure(pressure_pa)?;
        self.line.validate_jet_speed(speed)?;
        if self.card_centreline.is_some() {
            // Flux retains the windway's bulk speed, but convection, growth and Reynolds
            // validity use its centre speed. Preflight must use the proposed pressure,
            // including the developing exit's inviscid core, before any input is changed.
            let card_speed = if let Some(windway) = &self.unsteady_windway {
                windway.channel().steady_centre_speed_m_s(pressure_pa).map_err(FlueError::Physics)?
            } else { self.card_speed_at(pressure_pa, speed) };
            if (self.viscous_windway.is_some() || self.unsteady_windway.is_some())
                && card_speed > (2.0 * pressure_pa / self.rho).sqrt()
            {
                return Err(invalid("slit windway's steady centre speed exceeds the inviscid bound sqrt(2P/rho)"));
            }
            self.line.validate_jet_speed(card_speed)?;
        }
        Ok(speed)
    }
}

#[cfg(test)]
mod tests;

/// Scheduled physical outlet-pressure output through the shared audio interface.
pub mod performance;
