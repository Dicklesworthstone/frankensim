//! Supplied plate geometry, spatial closure, material memory and pressure phrase.
//! This is a bounded source adapter over the existing physical owners. Numeric
//! sections and supplied spectra remain authored data, not identified materials.
use super::{AperturePerformance, AperturePerformanceConfig, ApertureObservation, CoupledAperture};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use crate::bernoulli_aperture::dynamic::relaxation::{InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec};
use crate::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;
use crate::bernoulli_aperture::tube::{ApertureTube, UniformTubeSpec};
use crate::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::radiation::BaffledRadiationLoad;
use crate::render::schedule::GestureCompileError;
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;
use fs_material::{gas::GasState, visco::GeneralizedMaxwell};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use fs_scenario::gesture::GestureSchedule;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};

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
            Self::Geometry(e) => write!(f,"plate-valve geometry: {e}"),
            Self::Physics(e) => write!(f,"plate-valve physics: {e}"),
            Self::Gesture(e) => write!(f,"plate-valve gesture: {e}"),
        }
    }
}
impl std::error::Error for PlateValveInputError {}
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
    /// Actual source-node count.
    pub nodes: usize,
    /// Actual source-triangle count.
    pub triangles: usize,
    /// Distinct supplied section definitions, each used by the mesh.
    pub sections: usize,
    /// Positive retained material-memory coordinates.
    pub memory_branches: usize,
    /// Blowing-pressure input assignments compiled on the source clock.
    pub compiled_controls: usize,
    /// Length actually represented by the tube's integer transit [m].
    pub represented_tube_length_m: f64,
    /// Geometry-derived compact radiation replacing the memoryless terminal,
    /// absent on the original explicitly prescribed-reflection path.
    pub radiation_load: Option<BaffledRadiationLoad>,
}

/// Complete file-driven, geometry-derived valve and its finite scheduled runtime.
pub struct PlateValvePerformance {
    info: PlateValvePerformanceInfo,
    renderer: AperturePerformance,
}
impl PlateValvePerformance {
    /// Parse every record, derive the actual plate, bind its spatial lay and any
    /// supplied material memory, then compile the pressure phrase. Nothing steps
    /// during admission. Existing owners keep their physical/solver thresholds.
    /// The tube has lossless axial propagation. Its terminal is either an explicit
    /// memoryless reflection or a geometry-derived compact baffled radiation load.
    /// A selected radiation load replaces, never supplements, that reflection.
    /// Internal and exterior observations both use the accepted loaded mechanics.
    ///
    /// # Errors
    /// Invalid/count-limited source, cancellation, section/mesh/reduction/contact/
    /// material/propagation admission, or incomplete/unsupported gestures.
    pub fn from_bytes(bytes: &[u8], max_block: usize, gate: &CancelGate)
        -> Result<Self,PlateValveInputError>
    {
        checkpoint(gate)?;
        let parsed = Parsed::read(bytes,max_block,gate)?;
        let dt=1.0/f64::from(parsed.config.sample_rate_hz);
        // Reject an inadmissible radiation model before the plate eigensolve.
        // Both load and receiver use this same physical mouth and medium.
        let radiation_load=parsed.radiation_band.map(|band|BaffledRadiationLoad::new(
            parsed.tube.radius_m,parsed.air.density,parsed.air.sound_speed,dt,band,gate))
            .transpose().map_err(PlateValveInputError::Physics)?;
        let plate=PlateApertureReduction::from_chart(parsed.chart,parsed.plate,gate)
            .map_err(PlateValveInputError::Physics)?;
        if dt*(plate.stiffness_n_m()/plate.mass_kg()).sqrt()>parsed.max_angular_step {
            return Err(bad(0,"retained plate mode exceeds the declared mechanical angular-step allowance"));
        }
        let z=parsed.tube.characteristic_impedance(parsed.air.density).map_err(PlateValveInputError::Physics)?;
        let mut valve=DynamicAperture::from_plate_with_closure(plate,parsed.closure,parsed.air.density,z,
            dt,parsed.config.samples,parsed.initial).map_err(PlateValveInputError::Physics)?;
        if let Some((material,initial))=parsed.relaxation {
            valve=valve.with_plate_relaxation(material,initial).map_err(PlateValveInputError::Physics)?;
        }
        let memory_branches=valve.relaxation().map_or(0,|r|r.branches().len());
        let (system,observation,represented)=if let Some(load)=radiation_load {
            // A single existing network section supplies the SAME propagation
            // and physical inlet, with the existing stateful terminal owner.
            // No second bore, impedance solver or end-length correction is added.
            let t=parsed.tube;
            let network=ApertureNetwork::new(valve,TubeNetworkSpec {
                nodes:vec![NetworkNode::Inlet,load.termination()],
                sections:vec![TubeSection {nodes:[0,1],length_m:t.length_m,
                    radius_m:t.radius_m,max_length_error_m:t.max_length_error_m}],
                sound_speed_m_s:t.sound_speed_m_s,max_wave_memory_bytes:t.max_wave_memory_bytes,
            }).map_err(PlateValveInputError::Physics)?;
            let represented=network.represented_sections()[0].represented_length_m;
            let observation=match parsed.observation {
                ApertureObservation::Inlet=>ApertureObservation::Inlet,
                ApertureObservation::TubeTerminal=>ApertureObservation::NetworkNode(1),
                ApertureObservation::TubeBaffled(receiver)=>ApertureObservation::NetworkBaffled {node:1,receiver},
                _=>unreachable!("source grammar admits one inlet, terminal or outlet receiver"),
            };
            (CoupledAperture::Network(network),observation,represented)
        } else {
            // Preserve the original arithmetic and allocations when unselected.
            let tube=ApertureTube::new(valve,parsed.tube).map_err(PlateValveInputError::Physics)?;
            let represented=tube.represented_length_m();
            (CoupledAperture::Tube(tube),parsed.observation,represented)
        };
        let renderer=AperturePerformance::new(system,observation,parsed.schedule,parsed.config)
            .map_err(PlateValveInputError::Gesture)?;
        checkpoint(gate)?;
        let info=PlateValvePerformanceInfo {
            input_hash:hash_domain("org.frankensim.fs-couple.plate-valve-performance.v1",bytes),
            sample_rate_hz:parsed.config.sample_rate_hz,samples:parsed.config.samples,full_scale_pa:parsed.full_scale_pa,
            nodes:parsed.nodes,triangles:parsed.triangles,sections:parsed.sections,memory_branches,
            compiled_controls:renderer.pending_controls().len(),represented_tube_length_m:represented,radiation_load,
        };
        Ok(Self {info,renderer})
    }
    /// Exact source facts, retained separately from the evolving physical state.
    #[must_use]
    pub const fn info(&self) -> PlateValvePerformanceInfo {self.info}
    /// Complete finite renderer with immutable access to the physical specimen.
    #[must_use]
    pub const fn renderer(&self) -> &AperturePerformance {&self.renderer}
    /// Move the complete finite runtime into the existing pressure/PCM consumers.
    #[must_use]
    pub fn into_renderer(self) -> AperturePerformance {self.renderer}
}

struct Parsed {
    config:AperturePerformanceConfig,full_scale_pa:f64,air:GasState,tube:UniformTubeSpec,
    radiation_band:Option<f64>,observation:ApertureObservation,chart:PlateChart,plate:PlateApertureOptions,initial:ApertureState,
    closure:PlateClosureSpec,relaxation:Option<(PlateRelaxationSpec,InitialApertureMemory)>,
    max_angular_step:f64,schedule:GestureSchedule,nodes:usize,triangles:usize,sections:usize,
}
impl Parsed {
    fn read(bytes:&[u8],max_block:usize,gate:&CancelGate)->Result<Self,PlateValveInputError> {
        if bytes.len()>MAX_PLATE_VALVE_PERFORMANCE_BYTES || !(1..=65536).contains(&max_block) {
            return Err(bad(0,"input exceeds 4 MiB or callback capacity is outside 1..=65536"));
        }
        let text=std::str::from_utf8(bytes).map_err(|_|bad(0,"input must be UTF-8"))?;
        let (header,gestures)=text.split_once("\nschedule\n").ok_or_else(||bad(0,"missing LF-delimited schedule and canonical gesture suffix"))?;
        let mut r=Reader{lines:header.lines(),line:0};
        r.row(PLATE_VALVE_PERFORMANCE_SCHEMA)?.finish()?;
        let mut row=r.row("audio")?;
        let rate:u32=row.parse()?;let samples:u64=row.parse()?;let full_scale_pa=row.scalar()?;row.finish()?;
        if rate==0 || rate>192000 || samples==0 || samples>28800000 || samples>600*u64::from(rate) || full_scale_pa<=0.0 {
            return Err(bad(r.line,"audio requires positive rate <=192000, full scale and at most 600 seconds/28800000 samples"));
        }
        let mut row=r.row("ambient")?;
        let air=GasState::try_new_moist_air(row.scalar()?,row.scalar()?,row.scalar()?)
            .map_err(|_|bad(r.line,"ambient outside the shared moist-air domain"))?;row.finish()?;
        let mut row=r.row("tube")?;
        let length_m=row.scalar()?;let radius_m=row.scalar()?;
        let terminal=row.word()?;
        let (terminal_reflection,radiation_band)=if terminal=="baffled-low-ka" {
            // Zero here is only a placeholder in the geometry helper. It never
            // becomes a physical memoryless load: this selection uses a network.
            (0.0,Some(row.scalar()?))
        } else {
            let reflection:f64=terminal.parse().map_err(|_|bad(r.line,"tube requires a numeric reflectance or baffled-low-ka BAND_HZ"))?;
            if !reflection.is_finite() || reflection.abs()>1.0 {
                return Err(bad(r.line,"tube reflectance must be finite and passive"));
            }
            (reflection,None)
        };
        let tube=UniformTubeSpec{length_m,radius_m,terminal_reflection,
            max_length_error_m:row.scalar()?,max_wave_memory_bytes:row.count(64*1024*1024)?,sound_speed_m_s:air.sound_speed};row.finish()?;
        let mut row=r.row("observation")?;
        let observation=match row.word()? {
            "inlet"=>ApertureObservation::Inlet,"terminal"=>ApertureObservation::TubeTerminal,
            "baffled-outlet"=>ApertureObservation::TubeBaffled(CircularOutletReceiver {
                position_m:[row.scalar()?,row.scalar()?,row.scalar()?],
                radial_rings:row.count(128)?,angular_points:row.count(512)?,maximum_frequency_hz:row.scalar()?,
            }),
            _=>return Err(bad(r.line,"observation must explicitly name inlet, terminal or baffled-outlet")),
        };row.finish()?;
        if let (Some(band),ApertureObservation::TubeBaffled(receiver))=(radiation_band,observation) {
            if receiver.maximum_frequency_hz>band {
                return Err(bad(r.line,"exterior receiver band cannot exceed the selected radiation-load band"));
            }
        }
        let mut row=r.row("plate")?;
        let support=match row.word()? {"clamped"=>EdgeSupport::Clamped,"simply-supported"=>EdgeSupport::SimplySupported,
            _=>return Err(bad(r.line,"unsupported plate support law"))};
        let pretension=row.scalar()?;let damping_ratio=row.scalar()?;let low=row.scalar()?;let high=row.scalar()?;
        let mode_index=row.count(MAX_NODES*3)?;let max_angular_step=row.scalar()?;row.finish()?;
        if low<0.0 || high<=low || max_angular_step<=0.0 || max_angular_step>1.0 {
            return Err(bad(r.line,"plate requires an ordered nonnegative frequency window and angular step in (0,1]"));
        }
        let mut row=r.row("aperture")?;
        let rest_opening_m=row.scalar()?;let max_slit_mode_variation=row.scalar()?;let max_slope=row.scalar()?;row.finish()?;
        let mut row=r.row("initial")?;
        let initial=ApertureState{opening_m:row.scalar()?,opening_velocity_m_s:row.scalar()?};row.finish()?;
        let mut row=r.row("contact")?;
        let stiffness_pa_per_m_alpha=row.scalar()?;let alpha=row.scalar()?;
        let internal_loss_s_per_m=row.scalar()?;let max_penetration_m=row.scalar()?;row.finish()?;
        let mut row=r.row("source")?;let label=row.word()?.to_string();row.finish()?;
        if label.len()>256 {return Err(bad(r.line,"source label exceeds 256 bytes"));}
        let provenance=format!("supplied numeric plate-valve source {label}; no inferred material identification");
        let section_count=r.count("sections",MAX_SECTIONS)?;
        if section_count==0 {return Err(bad(r.line,"at least one material section is required"));}
        let mut sections=Vec::with_capacity(section_count);
        for _ in 0..section_count {
            let mut row=r.row("section")?;let kind=row.word()?;
            let thickness=row.scalar()?;let density=row.scalar()?;
            let section=match kind {
                "isotropic"=>PlateSection::isotropic(row.scalar()?,row.scalar()?,thickness,density),
                "orthotropic"=>PlateSection::orthotropic_plane_stress_at_angle(row.scalar()?,row.scalar()?,
                    row.scalar()?,row.scalar()?,thickness,density,row.scalar()?),
                _=>return Err(bad(r.line,"section must be explicitly isotropic or orthotropic")),
            }.map_err(PlateValveInputError::Geometry)?;row.finish()?;sections.push(section);
        }
        let node_count=r.count("nodes",MAX_NODES)?;
        if node_count<3 {return Err(bad(r.line,"at least three original mesh nodes are required"));}
        let mut nodes=Vec::with_capacity(node_count);let mut gaps=Vec::with_capacity(node_count);let mut supports=Vec::new();
        for i in 0..node_count {
            checkpoint(gate)?;let mut row=r.row("node")?;
            nodes.push((row.scalar()?,row.scalar()?));gaps.push(row.scalar()?);
            if row.flag()? {supports.push(i);}row.finish()?;
        }
        let triangle_count=r.count("triangles",MAX_TRIANGLES)?;
        if triangle_count==0 {return Err(bad(r.line,"at least one original mesh triangle is required"));}
        let mut triangles=Vec::with_capacity(triangle_count);let mut assigned=Vec::with_capacity(triangle_count);
        let mut section_triangles=vec![Vec::new();section_count];let mut lay_triangles=Vec::new();
        for i in 0..triangle_count {
            checkpoint(gate)?;let mut row=r.row("triangle")?;
            triangles.push([row.parse()?,row.parse()?,row.parse()?]);let section:usize=row.parse()?;
            assigned.push(*sections.get(section).ok_or_else(||bad(r.line,"triangle section does not exist"))?);
            section_triangles[section].push(i);if row.flag()? {lay_triangles.push(i);}row.finish()?;
        }
        if section_triangles.iter().any(Vec::is_empty) {return Err(bad(r.line,"every supplied section must be used by the plate"));}
        let edge_count=r.count("slit_edges",node_count)?;let mut slit_edges=Vec::with_capacity(edge_count);
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
        let mut row=r.row("compile_limits")?;
        let max_compile_work:u64=row.parse()?;let max_controls=row.count(MAX_CONTROLS)?;row.finish()?;
        if max_compile_work>16777216 {return Err(bad(r.line,"pressure compilation exceeds 16777216 visits"));}
        if r.lines.next().is_some() {return Err(bad(r.line+1,"unexpected trailing geometry record"));}
        let config=AperturePerformanceConfig{sample_rate_hz:rate,samples,max_block,max_compile_work,max_controls};
        let schedule=decode_schedule(gestures)?;
        Ok(Self{config,full_scale_pa,air,tube,radiation_band,observation,chart,plate,initial,closure,relaxation,max_angular_step,
            schedule,nodes:node_count,triangles:triangle_count,sections:section_count})
    }
}
fn decode_schedule(text:&str)->Result<GestureSchedule,PlateValveInputError> {
    let lines=text.lines().count();
    for line in text.lines() {
        if let Some(n)=line.strip_prefix("tracks\t") {
            if n.parse::<usize>().ok()!=Some(1) {return Err(bad(0,"one canonical pressure track is required"));}
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
    fn row(&mut self,key:&str)->Result<Row<'a>,PlateValveInputError> {
        self.line+=1;let line=self.lines.next().ok_or_else(||bad(self.line,"missing record"))?;
        let mut fields=line.split_ascii_whitespace();
        if fields.next()!=Some(key) {return Err(bad(self.line,"unexpected record kind or order"));}
        Ok(Row{fields,line:self.line})
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
}
