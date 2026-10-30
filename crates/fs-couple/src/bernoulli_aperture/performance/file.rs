//! Supplied plate geometry, spatial closure, material memory and pressure phrase.
//! This is a bounded source adapter over the existing physical owners. Numeric
//! sections and supplied spectra remain authored data, not identified materials.
//!
//! The aperture record accepts an optional fourth field: `single-blade` (the
//! legacy default) or `symmetric-opposed-pair`. For example,
//! `aperture 0.0002 0.25 0.1 symmetric-opposed-pair` selects two identical mirrored
//! blades without changing any solver. Mesh, sections and material history
//! describe one representative blade; openings and contact clearances remain
//! TOTAL gaps. Both blades use that same source and equal mirrored histories.
//! Force records, when present, prescribe TOTAL force equally split between
//! the blades. Asymmetric specimens/loads are outside this declared subspace.
//! An optional `slit_length LENGTH_M` record immediately after `aperture`
//! selects spatial viscous channel flow. Its viscosity comes from the actual
//! inlet gas, including regional composition and temperature declarations.
//!
//! A material section is either pack-backed,
//! `section material isotropic THICK PACK STATE [AXIS=VALUE ...]` or
//! `section material orthotropic12|orthotropic21 THICK ANGLE PACK STATE [AXIS=VALUE ...]`,
//! resolved through `crate::material_source` from a caller-supplied
//! material-card pack compiled at the named STATE; or bare numbers ending in
//! an explicit `fixture LABEL` (`section isotropic THICK RHO E NU fixture
//! LABEL`, `section orthotropic THICK RHO E1 E2 NU12 G12 ANGLE fixture
//! LABEL`). Bare numbers without a label refuse (bead pexfp.32). Each
//! section's pack, hash, state and query point, or its fixture label, is kept
//! in `section_provenance` for the render sidecar.
//!
//! After a DG source and before optional finite gas or observation records,
//! `dg_wall LAYOUT` selects the shared original material-wall grammar: rigid,
//! isotropic, orthotropic, physically positioned regions or a main/chimney
//! assembly. Its native solid state joins the same acoustic propagator before
//! the moving plate is bound. Wall and plate declarations share one immutable
//! artifact cache; their different material conditions remain independent.
//!
//! Optional `finite_wall solid` records follow finite gas, if present, and precede
//! observation. They select source-derived wall caloric properties, a finite
//! inventory and signed sample-addressed powers. Finite gas and wall together
//! require an explicit reciprocal conduction interface. The complete grammar
//! and its fixed-geometry all-solid scope are documented in `thermal_wall`.
mod force;
mod finite_gas;
mod thermal_wall;
mod wall;
/// Material-card multimode valves coupled to both ducts and one physical receiver.
pub mod material_modal;
use super::force::ApertureForceProgram;
use super::{AperturePerformance, AperturePerformanceConfig, ApertureObservation, CoupledAperture};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use crate::bernoulli_aperture::dynamic::relaxation::{InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec};
use crate::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction, PlateApertureSymmetry};
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;
use crate::bernoulli_aperture::tube::{ApertureTube, UniformTubeSpec};
use crate::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::radiation::BaffledRadiationLoad;
use crate::bernoulli_aperture::viscous::ViscousSlit;
use crate::render::schedule::GestureCompileError;
use crate::material_source::{MaterialRecord, fixture_label, plate_section};
use crate::thin_plate::PlateMaterialModel;
use fs_matdb::NormalizedMaterialCardPack;
use crate::dg_port::input::wall::regions::DgMaterialWallLayout;
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_blake3::{ContentHash, hash_domain};
use crate::bernoulli_aperture::viscothermal::ViscothermalSection;
use fs_exec::CancelGate;
use fs_material::{gas::GasState, visco::GeneralizedMaxwell};
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use fs_scenario::gesture::GestureSchedule;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
mod duct;
mod gas;
use duct::DuctInput;

/// Complete source grammar; the pressure suffix uses the existing gesture schema.
pub const PLATE_VALVE_PERFORMANCE_SCHEMA: &str = "frankensim-plate-valve-performance-v1";
/// Bound bytes before reading/decoding. Individual counts also have explicit caps.
pub const MAX_PLATE_VALVE_PERFORMANCE_BYTES: usize = 4 * 1024 * 1024;
const MAX_NODES: usize = 512;
const MAX_TRIANGLES: usize = 2048;
const MAX_SECTIONS: usize = 64;
const MAX_CONTROLS: usize = 262_144;

/// A refusal keeps its physical owner or one-based source line.
#[derive(Debug)]
pub enum PlateValveInputError {
    /// Malformed, unsupported or over-budget source record.
    Input { /// One-based header line, or zero for whole-file/gesture admission.
        line: usize, /// Failed admission condition.
        what: &'static str },
    /// Named gas admission, retaining the material owner and source line.
    Gas {
        /// One-based source line.
        line: usize,
        /// Original material/phase/composition refusal.
        source: fs_material::MaterialError,
    },
    /// Original stationary heat/species field refusal.
    SteadyGas { line: usize, source: fs_material::gas_field::steady::SteadyGasError },
    /// Original finite thermochemical inventory or caloric refusal.
    GasEnergy { line: usize, source: fs_material::thermochemical_gas::energy::GasEnergyError },
    /// Original source-derived wall caloric chart or enthalpy refusal.
    Phase { line: usize, source: fs_material::phase::PhaseStateError },
    /// Caller-supplied normalized material artifact could not be loaded.
    MaterialSource { /// One-based declaration line.
        line: usize, /// Original source refusal.
        source: String },
    /// Original source-bound mutual radiation artifact admission.
    RadiationSource { /// One-based performance declaration line.
        line: usize, /// Artifact loader, decoder or gas-owner refusal.
        source: crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError },
    /// Constitutive material-state resolution failed without a fallback.
    Material { /// One-based declaration line.
        line: usize, /// Original material owner refusal.
        source: fs_material::state_point::MaterialStatePointError },
    /// Mesh/section owner refusal.
    Geometry(fs_plate::PlateError),
    /// Reduction, coupled physics or material owner refusal.
    Physics(AcousticRealizeError),
    /// Existing pressure compiler refusal.
    Gesture(GestureCompileError),
}
impl core::fmt::Display for PlateValveInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f,"plate-valve input line {line}: {what}"),
            Self::Gas { line, source } => write!(f,"plate-valve gas line {line}: {source}"),
            Self::SteadyGas { line, source } => write!(f,"plate-valve steady gas line {line}: {source}"),
            Self::GasEnergy { line, source } => write!(f,"plate-valve finite gas line {line}: {source}"),
            Self::Phase { line, source } => write!(f,"plate-valve wall caloric line {line}: {source}"),
            Self::MaterialSource { line, source } => write!(f,"plate-valve material source line {line}: {source}"),
            Self::RadiationSource { line, source } => write!(f,"plate-valve radiation source line {line}: {source}"),
            Self::Material { line, source } => write!(f,"plate-valve material line {line}: {source}"),
            Self::Geometry(e) => write!(f,"plate-valve geometry: {e}"),
            Self::Physics(e) => write!(f,"plate-valve physics: {e}"),
            Self::Gesture(e) => write!(f,"plate-valve gesture: {e}"),
        }
    }
}
impl std::error::Error for PlateValveInputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Gas { source, .. } => Some(source),
            Self::SteadyGas { source, .. } => Some(source),
            Self::GasEnergy { source, .. } => Some(source),
            Self::Phase { source, .. } => Some(source),
            Self::Material { source, .. } => Some(source),
            Self::RadiationSource { source, .. } => Some(source),
            _ => None,
        }
    }
}
fn bad(line: usize, what: &'static str) -> PlateValveInputError { PlateValveInputError::Input { line, what } }
fn checkpoint(gate: &CancelGate) -> Result<(),PlateValveInputError> {
    if gate.is_requested() { return Err(bad(0,"plate-valve input cancelled")); } Ok(())
}

/// Source facts and physical reduction results; no material/fidelity authority.
#[derive(Clone, Copy, Debug)]
pub struct PlateValvePerformanceInfo {
    /// Domain-separated identity of all exact input bytes, independent of path.
    pub input_hash: ContentHash,
    /// Actual solver clock [Hz].
    pub sample_rate_hz: u32,
    /// Finite source duration in mechanical samples.
    pub samples: u64,
    /// PCM scale [Pa], not a per-part ensemble gain.
    pub full_scale_pa: f64,
    /// Actual source-node count on one representative blade, not duplicated nodes.
    pub nodes: usize,
    /// Actual source-triangle count on that representative blade.
    pub triangles: usize,
    /// Distinct supplied section definitions, each used by the mesh.
    pub sections: usize,
    /// Physical moving blades represented by the selected kinematic subspace.
    /// Two means identical mirrored blades, not independently controlled modes.
    pub moving_blades: usize,
    /// Positive retained material-memory coordinates.
    pub memory_branches: usize,
    /// Blowing-pressure input assignments compiled on the source clock.
    pub compiled_controls: usize,
    /// Explicit original-mesh force footprints.
    pub force_ports: usize,
    /// Raw sample-addressed mechanical force assignments, separate from pressure.
    pub force_controls: usize,
    /// Actual acoustic graph nodes, including the unique inlet and every load.
    pub duct_nodes: usize,
    /// Physical propagating sections; a legacy tube has one.
    pub duct_sections: usize,
    /// Count of geometry-derived radiating terminals, not receiver count.
    pub radiation_terminals: usize,
    /// Sum of represented physical section lengths [m], NOT a unique path
    /// length in a branched graph. Equal to the original tube value for one section.
    pub represented_tube_length_m: f64,
    /// Geometry-derived compact radiation replacing the memoryless terminal,
    /// absent on the original explicitly prescribed-reflection path. In a graph
    /// this describes the observed outlet (or a unique terminal for an internal
    /// observation). See `radiation_loads()` for every terminal.
    pub radiation_load: Option<BaffledRadiationLoad>,
}

/// Complete file-driven, geometry-derived valve and its finite scheduled runtime.
pub struct PlateValvePerformance {
    info: PlateValvePerformanceInfo,
    renderer: AperturePerformance,
    radiation_loads: Vec<(usize, BaffledRadiationLoad)>,
    viscothermal_losses: Vec<ViscothermalSection>,
    ambient_gas: GasState,
    exterior_gas: Option<GasState>,
    gas_records: Vec<GasRecord>,
    section_provenance: Vec<String>,
    dg_wall: Option<DgMaterialWallLayout>,
}
impl PlateValvePerformance {
    /// Parse every record, derive the actual plate and selected blade symmetry,
    /// bind its spatial contact and supplied material memory, then compile the
    /// pressure phrase. Nothing steps during admission. Existing owners keep
    /// their physical/solver thresholds. Omitting the aperture symmetry token
    /// preserves the original single-blade construction without a conversion.
    /// Legacy propagating intervals are lossless; an explicit viscothermal selection
    /// places gas-derived series/thermal memories along each selected section.
    /// `dg_tube` instead binds gas-derived Bessel memories inside the original DG
    /// exponential and solves its exact mean inlet relation with this same valve.
    /// Its computed baffled outlet also supports `observation baffled-outlet`,
    /// using the accepted terminal flow and that load's original exterior gas.
    /// Other graph loads retain explicit impedances, enclosed cavities, compliant walls
    /// and geometry-derived compact baffled radiation.
    /// A selected radiation load replaces, never supplements, that reflection.
    /// Internal and exterior observations both use the accepted loaded mechanics.
    ///
    /// # Errors
    /// Invalid/count-limited source, cancellation, section/mesh/reduction/contact/
    /// material/propagation admission, or incomplete/unsupported gestures.
    pub fn from_bytes(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
    ) -> Result<Self, PlateValveInputError> {
        Self::from_bytes_with_gas_sources(bytes, max_block, gate, None)
    }

    /// Decode every ambient, exterior and regional declaration against an
    /// explicit source registry, retaining the original gas records.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, PlateValveInputError> {
        Self::from_bytes_with_all_sources(bytes, max_block, gate, sources,
            |_| Err("a pack-backed material section needs an explicit material-card loader".into()))
    }

    /// As [`Self::from_bytes_with_gas_sources`], with `load_material` supplying
    /// each material-card pack a `section material` record names (by its path
    /// relative to the input). Each pack is loaded once.
    ///
    /// # Errors
    /// Those of [`Self::from_bytes`], a loader refusal, a card compiled at a
    /// state other than the record's, and a missing section property, by name.
    pub fn from_bytes_with_all_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlateValveInputError> {
        checkpoint(gate)?;
        let mut parsed = Parsed::read_with_gas_sources(bytes, max_block, gate, sources, &mut load_material)?;
        let dt=1.0/f64::from(parsed.config.sample_rate_hz);
        // One exterior declaration drives radiation; all internal media and
        // constitutive loads still come from their original regional states.
        let duct=parsed.duct.prepare_in_exterior_with_wall(&parsed.air,dt,parsed.observation,
            parsed.exterior.as_ref(),parsed.dg_wall.as_ref(),gate)?;
        let inlet_gas=duct.inlet_gas(parsed.air)?;
        let inlet_density=inlet_gas.density;
        let z=duct.inlet_impedance(inlet_density)?;
        let observation=duct.observation;
        let radiation_load=duct.observed_radiation();
        let radiation_loads=duct.radiation.clone();
        let viscothermal_losses=duct.losses.clone();
        let plate=PlateApertureReduction::from_chart(parsed.chart,parsed.plate,gate)
            .map_err(PlateValveInputError::Physics)?;
        let plate=match parsed.symmetry {
            PlateApertureSymmetry::SingleBlade => plate,
            PlateApertureSymmetry::SymmetricOpposedPair => plate.into_symmetric_opposed_pair()
                .map_err(PlateValveInputError::Physics)?,
        };
        let moving_blades=plate.moving_blade_count();
        if dt*(plate.stiffness_n_m()/plate.mass_kg()).sqrt()>parsed.max_angular_step {
            return Err(bad(0,"retained plate mode exceeds the declared mechanical angular-step allowance"));
        }
        // The opening coordinate and any declared Maxwell arms are propagated
        // exactly (DEF-04): this render path has no frequency-warping band.
        let mut valve=DynamicAperture::from_plate_with_closure_and_propagator(plate,parsed.closure,
            inlet_density,z,dt,parsed.config.samples,parsed.initial,
            crate::bernoulli_aperture::paired::BladePropagator::ExactZoh)
            .map_err(PlateValveInputError::Physics)?;
        if let Some((material,initial))=parsed.relaxation {
            valve=valve.with_plate_relaxation(material,initial).map_err(PlateValveInputError::Physics)?;
        }
        if let Some(length_m)=parsed.slit_length_m {
            let slit=ViscousSlit::from_gas(&inlet_gas,length_m)
                .map_err(PlateValveInputError::Physics)?;
            valve.set_viscous_slit(slit).map_err(PlateValveInputError::Physics)?;
        }
        let memory_branches=valve.relaxation().map_or(0,|r|r.branches().len());
        let system=duct.bind(valve)?;
        let (duct_nodes,duct_sections,represented)=match &system {
            CoupledAperture::Dg(d)=>(2,1,d.load().duct().total_length()),
            CoupledAperture::Tube(t)=>(2,1,t.represented_length_m()),
            CoupledAperture::Network(n)=>(n.spec().nodes.len(),n.spec().sections.len(),
                n.represented_sections().iter().map(|s|s.represented_length_m).sum()),
        };
        if !represented.is_finite() {return Err(bad(0,"total represented duct length overflowed"));}
        let dg_exterior = match &system {
            CoupledAperture::Dg(model) => model.load().radiation_record()
                .map(|record| *record.gas_record().state()),
            _ => None,
        };
        let dg_radiation_terminals = usize::from(dg_exterior.is_some());
        let mut renderer=match (parsed.exterior, observation) {
            (Some(exterior), ApertureObservation::NetworkBaffled { .. }) =>
                AperturePerformance::new_with_exterior_gas(system,observation,parsed.schedule,parsed.config,exterior),
            _ => AperturePerformance::new(system,observation,parsed.schedule,parsed.config),
        }.map_err(PlateValveInputError::Gesture)?;
        if let Some(program) = parsed.forces {
            renderer = renderer.with_plate_forces(program).map_err(PlateValveInputError::Gesture)?;
        }
        if let Some(wall) = parsed.finite_wall {
            renderer = wall.bind(renderer, parsed.finite_gas, gate, &mut parsed.gas_records)?;
        } else if let Some(heat) = parsed.finite_gas {
            renderer = heat.bind(renderer, gate, &mut parsed.gas_records)?;
        }
        checkpoint(gate)?;
        let info=PlateValvePerformanceInfo {
            input_hash:hash_domain("org.frankensim.fs-couple.plate-valve-performance.v1",bytes),
            sample_rate_hz:parsed.config.sample_rate_hz,samples:parsed.config.samples,full_scale_pa:parsed.full_scale_pa,
            nodes:parsed.nodes,triangles:parsed.triangles,sections:parsed.sections,moving_blades,memory_branches,
            force_ports:renderer.force_ports().len(),force_controls:renderer.pending_force_controls().len(),
            compiled_controls:renderer.pending_controls().len(),represented_tube_length_m:represented,radiation_load,
            duct_nodes,duct_sections,radiation_terminals:radiation_loads.len()+dg_radiation_terminals,
        };
        Ok(Self {
            info,
            renderer,
            radiation_loads,
            viscothermal_losses,
            ambient_gas: parsed.air,
            exterior_gas: dg_exterior.or(parsed.exterior),
            gas_records: parsed.gas_records,
            section_provenance: parsed.section_provenance,
            dg_wall: parsed.dg_wall,
        })
    }
    /// Each material section's source, in section order: `pack PATH HASH state
    /// STATE at AXIS=VALUE,...` for a pack-backed section, `fixture LABEL` for
    /// labelled bare numbers.
    #[must_use]
    pub fn section_provenance(&self) -> &[String] {&self.section_provenance}
    /// Exact source facts, retained separately from the evolving physical state.
    #[must_use]
    pub const fn info(&self) -> PlateValvePerformanceInfo {self.info}
    /// Resolved ambient/default material state. Explicit regional states remain
    /// available on the network; this is not an independent exterior medium.
    #[must_use]
    pub const fn ambient_gas(&self) -> &GasState {&self.ambient_gas}
    /// Explicit exterior used by every declared radiation load, including on
    /// internal-pressure observations. None preserves local terminal media.
    #[must_use]
    pub const fn exterior_gas(&self) -> Option<&GasState> {self.exterior_gas.as_ref()}
    /// Original gas declarations and source receipts, in authored read order.
    #[must_use]
    pub fn gas_records(&self) -> &[GasRecord] {
        &self.gas_records
    }
    /// Original authored DG material walls, including physical regions and
    /// explicit chimney assignments. None preserves the original omitted-wall path.
    #[must_use]
    pub const fn dg_wall_layout(&self) -> Option<&DgMaterialWallLayout> {self.dg_wall.as_ref()}
    /// Complete finite renderer with immutable access to the physical specimen.
    #[must_use]
    pub const fn renderer(&self) -> &AperturePerformance {&self.renderer}
    /// Compact radiation boundaries, paired with their source graph nodes.
    /// Computed DG radiation retains its original record on the DG load instead.
    /// Observing either kind adds no second reaction.
    #[must_use]
    pub fn radiation_loads(&self) -> &[(usize, BaffledRadiationLoad)] {&self.radiation_loads}
    /// Original section assignments and the exact retained RL/RC loss sites.
    /// Empty means no distributed gas loss was selected; it does not imply
    /// the absence of separately declared terminal, wall or cavity dissipation.
    #[must_use]
    pub fn viscothermal_losses(&self) -> &[ViscothermalSection] {&self.viscothermal_losses}
    /// Move the complete finite runtime into the existing pressure/PCM consumers.
    #[must_use]
    pub fn into_renderer(self) -> AperturePerformance {self.renderer}
}

struct Parsed {
    gas_records: Vec<GasRecord>,
    config: AperturePerformanceConfig,
    full_scale_pa: f64,
    air: GasState,
    duct: DuctInput,
    dg_wall: Option<DgMaterialWallLayout>,
    exterior: Option<GasState>,
    observation: ApertureObservation,
    chart: PlateChart,
    plate: PlateApertureOptions,
    symmetry: PlateApertureSymmetry,
    slit_length_m: Option<f64>,
    initial: ApertureState,
    closure: PlateClosureSpec,
    relaxation: Option<(PlateRelaxationSpec, InitialApertureMemory)>,
    forces: Option<ApertureForceProgram>,
    finite_gas: Option<finite_gas::FiniteGasInput>,
    finite_wall: Option<thermal_wall::FiniteWallInput>,
    max_angular_step: f64,
    schedule: GestureSchedule,
    nodes: usize,
    triangles: usize,
    sections: usize,
    section_provenance: Vec<String>,
}
impl Parsed {
    fn read(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
    ) -> Result<Self, PlateValveInputError> {
        Self::read_with_gas_sources(bytes, max_block, gate, None,
            &mut |_| Err("a pack-backed material section needs an explicit material-card loader".into()))
    }
    fn read_with_gas_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
        load_material: &mut dyn FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlateValveInputError> {
        if bytes.len()>MAX_PLATE_VALVE_PERFORMANCE_BYTES || !(1..=65536).contains(&max_block) {
            return Err(bad(0,"input exceeds 4 MiB or callback capacity is outside 1..=65536"));
        }
        let text=std::str::from_utf8(bytes).map_err(|_|bad(0,"input must be UTF-8"))?;
        let (header,gestures)=text.split_once("\nschedule\n").ok_or_else(||bad(0,"missing LF-delimited schedule and canonical gesture suffix"))?;
        let mut r=Reader{lines:header.lines(),line:0};
        r.row(PLATE_VALVE_PERFORMANCE_SCHEMA)?.finish()?;
        let mut row=r.row("audio")?;
        let rate:u32=row.parse()?;
        let samples:u64=row.parse()?;
        let full_scale_pa=row.scalar()?;
        row.finish()?;
        if rate==0 || rate>192000 || samples==0 || samples>28800000 || samples>600*u64::from(rate) || full_scale_pa<=0.0 {
            return Err(bad(r.line,"audio requires positive rate <=192000, full scale and at most 600 seconds/28800000 samples"));
        }
        let mut row=r.row("ambient")?;
        let mut gas_records = Vec::new();
        let air = gas::read_with_sources(&mut row, None, sources, &mut gas_records)?;
        row.finish()?;
        let exterior=if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("exterior") {
            let mut row=r.row("exterior")?;
            let state=gas::read_with_sources(&mut row,None,sources,&mut gas_records)?;
            if state.pressure != air.pressure {
                return Err(bad(row.line,"exterior and ambient must share static pressure; mean interface flow is not modeled"));
            }
            row.finish()?;Some(state)
        } else {None};
        let duct = DuctInput::read_with_gas_sources(&mut r, &air, gate, sources, &mut gas_records)?;
        let mut packs = wall::MaterialPacks::new();
        let dg_wall = wall::read(&mut r, &duct, &mut packs, load_material, gate)?;
        let finite_gas = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("finite_gas") {
            Some(finite_gas::FiniteGasInput::read(&mut r, sources, gate)?)
        } else { None };
        let finite_wall = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("finite_wall") {
            Some(thermal_wall::FiniteWallInput::read(&mut r, dg_wall.as_ref(), gate)?)
        } else { None };
        if finite_wall.as_ref().is_some_and(|wall| wall.requires_gas() != finite_gas.is_some()) {
            return Err(bad(r.line,
                "finite gas and wall heat require an explicit reciprocal conduction interface"));
        }
        let mut row=r.row("observation")?;
        let observation=match row.word()? {
            "inlet"=>ApertureObservation::Inlet,"terminal"=>ApertureObservation::TubeTerminal,
            "network-node"=>ApertureObservation::NetworkNode(row.parse()?),
            "baffled-outlet"=>ApertureObservation::TubeBaffled(read_receiver(&mut row)?),
            "network-baffled"=>ApertureObservation::NetworkBaffled {
                node:row.parse()?,receiver:read_receiver(&mut row)?,
            },
            _=>return Err(bad(r.line,"observation must name a physical inlet, terminal, network node or outlet receiver")),
        };
        row.finish()?;
        let mut row=r.row("plate")?;
        let support=match row.word()? {"clamped"=>EdgeSupport::Clamped,"simply-supported"=>EdgeSupport::SimplySupported,
            _=>return Err(bad(r.line,"unsupported plate support law"))};
        let pretension=row.scalar()?;
        let damping_ratio=row.scalar()?;
        let low=row.scalar()?;
        let high=row.scalar()?;
        let mode_index=row.count(MAX_NODES*3)?;
        let max_angular_step=row.scalar()?;
        row.finish()?;
        if low<0.0 || high<=low || max_angular_step<=0.0 || max_angular_step>1.0 {
            return Err(bad(r.line,"plate requires an ordered nonnegative frequency window and angular step in (0,1]"));
        }
        let mut row=r.row("aperture")?;
        let rest_opening_m=row.scalar()?;
        let max_slit_mode_variation=row.scalar()?;
        let max_slope=row.scalar()?;
        let symmetry=match row.fields.next() {
            None | Some("single-blade") => PlateApertureSymmetry::SingleBlade,
            Some("symmetric-opposed-pair") => PlateApertureSymmetry::SymmetricOpposedPair,
            Some(_) => return Err(bad(r.line,"aperture symmetry must be single-blade or symmetric-opposed-pair")),
        };
        row.finish()?;
        let slit_length_m=if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("slit_length") {
            let mut row=r.row("slit_length")?;
            let length_m=row.scalar()?;
            row.finish()?;
            if length_m<=0.0 {
                return Err(bad(r.line,"viscous slit channel length must be positive"));
            }
            Some(length_m)
        } else {None};
        let mut row=r.row("initial")?;
        let initial=ApertureState{opening_m:row.scalar()?,opening_velocity_m_s:row.scalar()?};
        row.finish()?;
        let mut row=r.row("contact")?;
        let stiffness_pa_per_m_alpha=row.scalar()?;
        let alpha=row.scalar()?;
        let internal_loss_s_per_m=row.scalar()?;
        let max_penetration_m=row.scalar()?;
        row.finish()?;
        let mut row=r.row("source")?;
        let label=row.word()?.to_string();
        row.finish()?;
        if label.len()>256 {return Err(bad(r.line,"source label exceeds 256 bytes"));}
        let provenance=format!("supplied numeric plate-valve source {label}; no inferred material identification");
        let section_count=r.count("sections",MAX_SECTIONS)?;
        if section_count==0 {return Err(bad(r.line,"at least one material section is required"));}
        let mut sections=Vec::with_capacity(section_count);
        let mut section_provenance=Vec::with_capacity(section_count);
        for _ in 0..section_count {
            checkpoint(gate)?;
            let mut row=r.row("section")?;let kind=row.word()?;let line=row.line;
            let ingress=|e:crate::material_source::MaterialIngressError|PlateValveInputError::MaterialSource{line,source:e.0};
            let (section,origin)=match kind {
                "isotropic"|"orthotropic"=>{
                    let thickness=row.scalar()?;let density=row.scalar()?;
                    let section=if kind=="isotropic" {PlateSection::isotropic(row.scalar()?,row.scalar()?,thickness,density)}
                        else {PlateSection::orthotropic_plane_stress_at_angle(row.scalar()?,row.scalar()?,
                            row.scalar()?,row.scalar()?,thickness,density,row.scalar()?)}.map_err(PlateValveInputError::Geometry)?;
                    let label=fixture_label(&row.rest()).map_err(ingress)?;
                    (section,format!("fixture {label}"))
                }
                "material"=>{
                    let model=match row.word()? {
                        "isotropic"=>PlateMaterialModel::Isotropic,
                        "orthotropic12"=>PlateMaterialModel::Orthotropic12,
                        "orthotropic21"=>PlateMaterialModel::Orthotropic21,
                        _=>return Err(bad(line,"a material section needs isotropic, orthotropic12 or orthotropic21")),
                    };
                    let thickness=row.scalar()?;
                    let angle=if model==PlateMaterialModel::Isotropic {0.0} else {row.scalar()?};
                    let record=MaterialRecord::parse(&row.rest()).map_err(ingress)?;
                    let pack=wall::load_pack(&mut packs,&record.pack,line,load_material,gate)?;
                    let card=record.card(&pack).map_err(ingress)?;
                    let point=record.query_point(card).map_err(ingress)?;
                    let (section,_)=plate_section(card,&point,model,thickness,angle).map_err(ingress)?;
                    (section,record.provenance(&pack))
                }
                _=>return Err(bad(line,"section must be isotropic, orthotropic or a pack-backed material section")),
            };
            sections.push(section);section_provenance.push(origin);
        }
        let node_count=r.count("nodes",MAX_NODES)?;
        if node_count<3 {return Err(bad(r.line,"at least three original mesh nodes are required"));}
        let mut nodes=Vec::with_capacity(node_count);
        let mut gaps=Vec::with_capacity(node_count);
        let mut supports=Vec::new();
        for i in 0..node_count {
            checkpoint(gate)?;let mut row=r.row("node")?;
            nodes.push((row.scalar()?,row.scalar()?));gaps.push(row.scalar()?);
            if row.flag()? {supports.push(i);}row.finish()?;
        }
        let triangle_count=r.count("triangles",MAX_TRIANGLES)?;
        if triangle_count==0 {return Err(bad(r.line,"at least one original mesh triangle is required"));}
        let mut triangles=Vec::with_capacity(triangle_count);
        let mut assigned=Vec::with_capacity(triangle_count);
        let mut section_triangles=vec![Vec::new();section_count];
        let mut lay_triangles=Vec::new();
        for i in 0..triangle_count {
            checkpoint(gate)?;let mut row=r.row("triangle")?;
            triangles.push([row.parse()?,row.parse()?,row.parse()?]);let section:usize=row.parse()?;
            assigned.push(*sections.get(section).ok_or_else(||bad(r.line,"triangle section does not exist"))?);
            section_triangles[section].push(i);if row.flag()? {lay_triangles.push(i);}row.finish()?;
        }
        if section_triangles.iter().any(Vec::is_empty) {return Err(bad(r.line,"every supplied section must be used by the plate"));}
        let edge_count=r.count("slit_edges",node_count)?;
        let mut slit_edges=Vec::with_capacity(edge_count);
        for _ in 0..edge_count {let mut row=r.row("edge")?;slit_edges.push([row.parse()?,row.parse()?]);row.finish()?;}
        let mesh=PlateMesh::from_unstructured(nodes,triangles).map_err(PlateValveInputError::Geometry)?;
        let chart=PlateChart::with_boundary_and_regions(mesh,sections[0],supports,Vec::new())
            .and_then(|c|c.with_element_sections(assigned)).map_err(PlateValveInputError::Geometry)?;
        let plate=PlateApertureOptions {assembly:AssemblyOptions{pretension,support},
            eigenvalue_window:((core::f64::consts::TAU*low).powi(2),(core::f64::consts::TAU*high).powi(2)),
            mode_index,eigensolver:SliceOptions::default(),max_nodes:MAX_NODES,max_triangles:MAX_TRIANGLES,
            slit_edges,rest_opening_m,damping_ratio,max_slit_mode_variation,max_slope};
        let closure=PlateClosureSpec{nodal_rest_gap_m:gaps,lay_triangles,stiffness_pa_per_m_alpha,alpha,
            internal_loss_s_per_m,provenance:provenance.clone(),max_penetration_m};
        let regions=r.count("relaxation",section_count)?;
        let relaxation=if regions==0 {None} else {
            if regions!=section_count {return Err(bad(r.line,"relaxation must cover every material section"));}
            let mut row=r.row("time_limits")?;let max_dt_over_tau=row.scalar()?;let memory_angle=row.scalar()?;row.finish()?;
            if memory_angle>max_angular_step {return Err(bad(r.line,"material angular-step allowance cannot exceed the declared mechanics allowance"));}
            let mut row=r.row("memory_initial")?;
            let initial=match row.word()? {
                "relaxed"=>InitialApertureMemory::Relaxed,"unrelaxed"=>InitialApertureMemory::Unrelaxed,
                "explicit"=>{let count=row.count(64)?;let mut values=Vec::with_capacity(count);
                    for _ in 0..count {values.push(row.scalar()?);}InitialApertureMemory::ViscousDisplacementM(values)},
                _=>return Err(bad(r.line,"memory must be explicitly relaxed, unrelaxed or supplied per arm")),
            };row.finish()?;
            let mut seen=vec![false;section_count];let mut maps=Vec::with_capacity(regions);let mut total=0;
            for _ in 0..regions {
                let mut row=r.row("region")?;let section:usize=row.parse()?;
                if section>=section_count || seen[section] {return Err(bad(r.line,"memory region must name a unique existing section"));}
                seen[section]=true;
                let e_inf=row.scalar()?;let poisson_ratio=row.scalar()?;let band_hz=(row.scalar()?,row.scalar()?);
                let branches=row.count(64-total)?;total+=branches;row.finish()?;
                let mut terms=Vec::with_capacity(branches);
                for _ in 0..branches {let mut row=r.row("branch")?;terms.push((row.scalar()?,row.scalar()?));row.finish()?;}
                let material=GeneralizedMaxwell::new(e_inf,terms).map_err(|_|bad(r.line,"invalid supplied Maxwell material"))?;
                maps.push(PlateRelaxationRegion{triangles:section_triangles[section].clone(),material,
                    poisson_ratio,band_hz,provenance:provenance.clone()});
            }
            Some((PlateRelaxationSpec{regions:maps,max_branches:64,max_dt_over_tau,max_angular_step:memory_angle},initial))
        };
        let mut row=r.next()?;
        let forces = match row.word()? {
            "compile_limits" => None,
            "force_ports" => {
                let count = row.count(32)?; row.finish()?;
                let program = force::read(&mut r,count,node_count,triangle_count,samples,gate)?;
                row = r.row("compile_limits")?;
                Some(program)
            }
            _ => return Err(bad(r.line,"expected force_ports or compile_limits")),
        };
        let max_compile_work:u64=row.parse()?;
        let max_controls=row.count(MAX_CONTROLS)?;
        row.finish()?;
        if max_compile_work>16777216 {return Err(bad(r.line,"pressure compilation exceeds 16777216 visits"));}
        if forces.as_ref().is_some_and(|f| f.events.len() > max_controls) {
            return Err(bad(r.line,"force assignments exceed the declared control budget"));
        }
        if r.lines.next().is_some() {return Err(bad(r.line+1,"unexpected trailing geometry record"));}
        let config=AperturePerformanceConfig{sample_rate_hz:rate,samples,max_block,max_compile_work,max_controls};
        let schedule=decode_schedule(gestures)?;
        Ok(Self{gas_records,config,full_scale_pa,air,duct,dg_wall,exterior,observation,chart,plate,symmetry,slit_length_m,initial,closure,relaxation,forces,finite_gas,finite_wall,max_angular_step,
            schedule,nodes:node_count,triangles:triangle_count,sections:section_count,section_provenance})
    }
}
fn read_receiver(row:&mut Row<'_>)->Result<CircularOutletReceiver,PlateValveInputError> {
    Ok(CircularOutletReceiver {position_m:[row.scalar()?,row.scalar()?,row.scalar()?],
        radial_rings:row.count(128)?,angular_points:row.count(512)?,maximum_frequency_hz:row.scalar()?})
}
fn decode_schedule(text:&str)->Result<GestureSchedule,PlateValveInputError> {
    decode_schedule_with_track_limit(text, 1)
}
fn decode_schedule_with_track_limit(text:&str,max_tracks:usize)->Result<GestureSchedule,PlateValveInputError> {
    let lines=text.lines().count();
    for line in text.lines() {
        if let Some(n)=line.strip_prefix("tracks\t") {
            if n.parse::<usize>().ok().is_none_or(|count| count==0 || count>max_tracks) {
                return Err(bad(0,"canonical tracks exceed the declared source consumer bound"));
            }
        }
        if let Some(n)=line.strip_prefix("events\t") {
            let count=n.parse::<usize>().map_err(|_|bad(0,"invalid gesture event count"))?;
            if count>16384 || count>lines {return Err(bad(0,"gesture event count exceeds source/16384-event budget"));}
        }
    }
    let schedule=GestureSchedule::from_canonical_bytes(text.as_bytes())
        .map_err(|e|PlateValveInputError::Gesture(GestureCompileError::Gesture(e)))?;
    if schedule.to_canonical_bytes()!=text.as_bytes() {return Err(bad(0,"gesture bytes must be canonical without ignored suffixes"));}
    Ok(schedule)
}
struct Reader<'a>{lines:Lines<'a>,line:usize}
struct Row<'a>{fields:SplitAsciiWhitespace<'a>,line:usize}
impl<'a> Reader<'a> {
    fn next(&mut self)->Result<Row<'a>,PlateValveInputError> {
        self.line+=1;let line=self.lines.next().ok_or_else(||bad(self.line,"missing record"))?;
        Ok(Row{fields:line.split_ascii_whitespace(),line:self.line})
    }
    fn row(&mut self,key:&str)->Result<Row<'a>,PlateValveInputError> {
        let mut row=self.next()?;
        if row.word()?!=key {return Err(bad(self.line,"unexpected record kind or order"));}Ok(row)
    }
    fn count(&mut self,key:&str,max:usize)->Result<usize,PlateValveInputError> {
        let mut row=self.row(key)?;let count=row.count(max)?;row.finish()?;Ok(count)
    }
}
impl<'a> Row<'a> {
    fn word(&mut self)->Result<&'a str,PlateValveInputError> {self.fields.next().ok_or_else(||bad(self.line,"missing field"))}
    fn parse<T:FromStr>(&mut self)->Result<T,PlateValveInputError> {self.word()?.parse().map_err(|_|bad(self.line,"invalid numeric field"))}
    fn scalar(&mut self)->Result<f64,PlateValveInputError> {
        let x:f64=self.parse()?;if !x.is_finite() {return Err(bad(self.line,"physical scalars must be finite"));}Ok(x)
    }
    fn count(&mut self,max:usize)->Result<usize,PlateValveInputError> {
        let n=self.parse()?;if n>max {return Err(bad(self.line,"count exceeds the input budget"));}Ok(n)
    }
    fn flag(&mut self)->Result<bool,PlateValveInputError> {
        match self.word()? {"0"=>Ok(false),"1"=>Ok(true),_=>Err(bad(self.line,"flag must be exactly 0 or 1"))}
    }
    fn finish(mut self)->Result<(),PlateValveInputError> {
        if self.fields.next().is_some() {return Err(bad(self.line,"unexpected extra field"));}Ok(())
    }
    /// Every remaining field, consuming the row.
    fn rest(&mut self)->Vec<&'a str> {self.fields.by_ref().collect()}
}

#[cfg(test)]
mod opposed_performance_tests {
    use super::*;
    use crate::pcm_wav::observation::PressureRenderer;

    const INPUT: &str = include_str!("../../../examples/plate-valve.performance");
    const APERTURE: &str = "aperture 0.0002 0.25 0.1";

    fn selected(selector: &str) -> String {
        assert_eq!(INPUT.matches(APERTURE).count(), 1);
        INPUT.replace(APERTURE, &format!("{APERTURE} {selector}"))
    }
    fn load(text: &str, block: usize) -> PlateValvePerformance {
        let resolved = text.replace("audio 96000 9602", "audio 192000 19204");
        PlateValvePerformance::from_bytes(resolved.as_bytes(), block, &CancelGate::new()).unwrap()
    }
    fn render(text: &str, block: usize) -> Vec<f64> {
        let voice = load(text, block);
        let mut out = vec![0.0; voice.info().samples as usize];
        let mut renderer = voice.into_renderer();
        for chunk in out.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        out
    }
    fn same_bits(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) { assert_eq!(x.to_bits(), y.to_bits()); }
    }
    fn different_pressure(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        assert!(a.iter().chain(b).all(|v| v.is_finite()));
        assert!(a.iter().zip(b).any(|(x, y)| (x-y).abs() > 1e-8));
    }

    #[test]
    fn single_blade_selector_preserves_legacy_pressure_bits() {
        let explicit = selected("single-blade");
        assert_eq!(load(INPUT, 37).info().moving_blades, 1);
        assert_eq!(load(&explicit, 64).info().moving_blades, 1);
        same_bits(&render(INPUT, 37), &render(&explicit, 64));
    }

    #[test]
    fn opposed_source_selection_changes_loaded_physics_and_preserves_callback_partitioning() {
        let text = selected("symmetric-opposed-pair");
        let single = load(INPUT, 37);
        let pair = load(&text, 113);
        assert_eq!(pair.info().moving_blades, 2);
        assert_eq!((pair.info().nodes, pair.info().triangles, pair.info().sections), (21, 24, 1));
        assert_eq!(pair.info().memory_branches, single.info().memory_branches);
        let a = single.renderer().system().aperture();
        let b = pair.renderer().system().aperture();
        assert_eq!(b.spec().mass_kg, 0.5*a.spec().mass_kg);
        assert_eq!(b.spec().stiffness_n_m, 0.5*a.spec().stiffness_n_m);
        assert_eq!(b.spec().aperture.closing_pressure_pa, 0.5*a.spec().aperture.closing_pressure_pa);
        assert_eq!(b.spec().aperture.rest_opening_m, a.spec().aperture.rest_opening_m);
        assert_eq!(b.spec().aperture.width_m, a.spec().aperture.width_m);
        assert_eq!(b.spec().density_kg_m3, a.spec().density_kg_m3);
        assert_eq!(b.spec().impedance_pa_s_m3, a.spec().impedance_pa_s_m3);
        assert_eq!(b.state(), a.state());
        let source = b.plate_reduction().unwrap();
        assert_eq!(source.symmetry(), PlateApertureSymmetry::SymmetricOpposedPair);
        assert!(source.material_chart().is_none(), "numeric sections cannot mint material receipts");
        assert_eq!(source.chart().mesh.nodes, a.plate_reduction().unwrap().chart().mesh.nodes);
        for (one, two) in a.relaxation().unwrap().branches().iter().zip(b.relaxation().unwrap().branches()) {
            assert!((two.stiffness/(0.5*one.stiffness)-1.0).abs() < 1e-12);
            assert_eq!(two.relaxation_time_s, one.relaxation_time_s);
        }
        let pressure = render(&text, 37);
        same_bits(&pressure, &render(&text, 113));
        different_pressure(&render(INPUT, 37), &pressure);
    }

    #[test]
    fn opposed_material_density_and_ambient_pressure_changes_reach_the_same_pressure_renderer() {
        let text = selected("symmetric-opposed-pair");
        let heavy = text.replace("0.0003 900 4000000000", "0.0003 1800 4000000000");
        let thinner_air = text.replace("ambient 293.15 101325 0", "ambient 293.15 90000 0");
        assert_ne!(heavy, text);
        assert_ne!(thinner_air, text);
        let a = load(&text, 37);
        let b = load(&heavy, 37);
        let c = load(&thinner_air, 37);
        let original = a.renderer().system().aperture().spec();
        let changed_mass = b.renderer().system().aperture().spec();
        let changed_gas = c.renderer().system().aperture().spec();
        assert!((changed_mass.mass_kg/original.mass_kg-2.0).abs() < 2e-5);
        assert_eq!(changed_gas.mass_kg, original.mass_kg);
        assert!(changed_gas.density_kg_m3 < original.density_kg_m3);
        assert_eq!(changed_gas.aperture.rest_opening_m, original.aperture.rest_opening_m);
        let pressure = render(&text, 37);
        different_pressure(&pressure, &render(&heavy, 37));
        different_pressure(&pressure, &render(&thinner_air, 37));
    }

    #[test]
    fn unsupported_or_repeated_blade_selectors_are_not_silently_downgraded() {
        for selector in ["independent-pair", "unknown", "2", "single-blade extra",
            "symmetric-opposed-pair extra", "symmetric-opposed-pair symmetric-opposed-pair"]
        {
            assert!(PlateValvePerformance::from_bytes(selected(selector).as_bytes(), 37,
                &CancelGate::new()).is_err(), "accepted {selector}");
        }
        let gate = CancelGate::new();
        gate.request();
        assert!(PlateValvePerformance::from_bytes(selected("symmetric-opposed-pair").as_bytes(),
            37, &gate).is_err());
    }
}
