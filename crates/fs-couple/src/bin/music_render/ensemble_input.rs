//! File-driven multi-rate physical ensembles; use existing model and PCM owners.
use std::io::{Read, Write};
use std::path::PathBuf;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::gas_species::input::ThermochemicalSources;
use fs_couple::bernoulli_aperture::performance::file::{PlateValvePerformance, PLATE_VALVE_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::performance::file::material_modal::{
    MaterialModalValvePerformance, MATERIAL_MODAL_VALVE_SCHEMA,
};
use fs_couple::bernoulli_aperture::performance::file::material_modal::dg::{
    MaterialModalDgValvePerformance, MATERIAL_MODAL_DG_VALVE_SCHEMA,
};
use fs_couple::bernoulli_aperture::radiation::environment::multiport::input::MAX_MIMO_RADIATION_BYTES;
use fs_couple::bernoulli_aperture::two_sided::ValveSide;
use fs_couple::bowed_string::runtime::schedule::file::{BowedPerformance, BOWED_PERFORMANCE_SCHEMA, MAX_BOWED_PERFORMANCE_BYTES};
use fs_couple::brass_loop::performance::{LipDescription, LIP_PERFORMANCE_SCHEMA, MAX_LIP_PERFORMANCE_BYTES};
use fs_couple::voice_performance::{VoiceDescription, VOICE_PERFORMANCE_SCHEMA, MAX_VOICE_PERFORMANCE_BYTES};
use fs_couple::flue_loop::dg::performance::file::{FlueDgFilePerformance, FLUE_DG_PERFORMANCE_SCHEMA};
use fs_couple::flue_loop::network::performance::file::{
    FlueFilePerformance, FLUE_PERFORMANCE_SCHEMA, MAX_FLUE_PERFORMANCE_BYTES,
};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::observation::ensemble::{PressureEnsemble, PressureEnsembleConfig};
use fs_couple::render::plate::file::{PlatePerformance, PLATE_PERFORMANCE_SCHEMA};
use fs_couple::render::schedule::{ScheduledRenderer, force::file::ModalPerformance};
use fs_couple::render::schedule::reed::{ReedPerformance, REED_PERFORMANCE_SCHEMA, MAX_REED_PERFORMANCE_BYTES};
use super::{RATE, create_outputs, json_string, stream_output};

// These aggregate limits supplement, never replace, each source loader's caps.
const MAX_PARTS: usize = 16;
const MAX_PART_BYTES: usize = 4 * 1024 * 1024;
const MAX_COMPONENTS: usize = 128;
const MAX_MODES: usize = 8192;
const MAX_CONTROLS: usize = 262_144;
const MAX_SAMPLES: u64 = 600 * RATE as u64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind { Modal, Plate, Bow, Reed, Valve, Flue, Lip, Voice }
impl Kind {
    fn label(self) -> &'static str {
        match self { Self::Modal => "modal", Self::Plate => "plate", Self::Bow => "bow", Self::Reed => "reed", Self::Valve => "valve", Self::Flue => "flue", Self::Lip => "lip", Self::Voice => "voice" }
    }
}
struct Input { kind: Kind, path: PathBuf }
struct Options {
    output: PathBuf,
    inputs: Vec<Input>,
    block: usize,
    full_scale_pa: f64,
    decimate: bool,
}
fn options(args: &[String]) -> Result<Options, String> {
    let mut output = None;
    let mut inputs = Vec::new();
    let mut block = None;
    let mut full_scale_pa = None;
    let mut decimate = false;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--modal" | "--plate" | "--bow" | "--reed" | "--valve" | "--flue" | "--lip" | "--voice" => {
                if inputs.len() == MAX_PARTS { return Err("ensemble exceeds the 16-part budget".into()); }
                let path = iter.next().filter(|s| !s.starts_with("--"))
                    .ok_or_else(|| format!("{arg} requires a source performance path"))?;
                inputs.push(Input { kind: match arg.as_str() { "--modal" => Kind::Modal, "--plate" => Kind::Plate, "--bow" => Kind::Bow, "--reed" => Kind::Reed, "--valve" => Kind::Valve, "--flue" => Kind::Flue, "--lip" => Kind::Lip, _ => Kind::Voice },
                    path: PathBuf::from(path) });
            }
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                block = Some(iter.next().and_then(|s| s.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block requires an integer in 1..=65536".to_string())?);
            }
            "--full-scale-pa" => {
                if full_scale_pa.is_some() { return Err("--full-scale-pa may only be supplied once".into()); }
                full_scale_pa = Some(iter.next().and_then(|s| s.parse::<f64>().ok())
                    .filter(|v| v.is_finite() && *v > 0.0)
                    .ok_or_else(|| "--full-scale-pa requires finite positive pascals".to_string())?);
            }
            "--decimate" => {
                if decimate { return Err("--decimate may only be supplied once".into()); }
                decimate = true;
            }
            value if value.starts_with('-') => return Err(format!("unsupported ensemble option {value:?}")),
            value => {
                if output.is_some() { return Err("ensemble accepts exactly one output path; inputs require --modal/--plate/--bow/--reed/--valve/--flue/--lip/--voice".into()); }
                output = Some(PathBuf::from(value));
            }
        }
    }
    if inputs.is_empty() { return Err("ensemble requires at least one explicitly typed source performance".into()); }
    Ok(Options {
        output: output.ok_or_else(|| "ensemble requires an output WAV path".to_string())?,
        inputs, block: block.unwrap_or(512), decimate,
        full_scale_pa: full_scale_pa.ok_or_else(|| "ensemble requires its own explicit --full-scale-pa; source PCM scales are not gains".to_string())?,
    })
}

struct SourceInfo {
    kind: Kind,
    schema: &'static str,
    hash: ContentHash,
    // Actual pressure-producer clock; an existing source may already own a
    // multirate microphone. The physical mechanics are reported separately.
    rate: u32,
    samples: u64,
    mechanics_rate: u32,
    mechanics_samples: u64,
    // Existing receiver's numerical latency on the producer clock. This never
    // includes its physical source-to-listener travel time.
    source_filter_delay_samples: f64,
    source_full_scale_pa: f64,
    components: usize,
    modes: usize,
    controls: usize,
    exterior_json: String,
    listener: Listener,
}
/// The gas a part's pressure is heard in, read from its own declaration.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Listener {
    Gas(ListenerMedium),
    /// The part's sound depends on no exterior gas; the reason is receipted.
    Unconstrained(&'static str),
}
#[derive(Clone, Copy, Debug, PartialEq)]
struct ListenerMedium { density: f64, sound_speed: f64, gas: Option<GasState> }
impl ListenerMedium {
    fn of(gas: &GasState) -> Self { Self { density: gas.density, sound_speed: gas.sound_speed, gas: Some(*gas) } }
    fn provenance_json(&self) -> String {
        match self.gas {
            Some(gas) => format!(
                "{{\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"state_authority\":\"complete_gas_state\",\"temperature_k\":{:e},\"pressure_pa\":{:e},\"dynamic_viscosity_pa_s\":{:e},\"thermal_conductivity_w_m_k\":{:e},\"heat_capacity_ratio\":{:e},\"specific_gas_constant_j_kg_k\":{:e},\"specific_heat_cp_j_kg_k\":{:e},\"prandtl\":{:e},\"characteristic_impedance_pa_s_m\":{:e},\"water_mole_fraction\":{:e}}}",
                gas.density, gas.sound_speed, gas.temperature, gas.pressure,
                gas.dynamic_viscosity, gas.thermal_conductivity, gas.gamma,
                gas.specific_gas_constant, gas.specific_heat_cp, gas.prandtl,
                gas.characteristic_impedance, gas.water_mole_fraction,
            ),
            None => format!(
                "{{\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"state_authority\":\"receiver_coefficients_only\"}}",
                self.density, self.sound_speed,
            ),
        }
    }
    /// The observer consumes rho (amplitude) and c (time of flight) bit for bit.
    /// If either source claims a complete gas state, both must retain and match
    /// that complete state. A receiver-only (rho, c) tuple cannot prove that its
    /// composition, temperature and pressure are the same as a sourced medium.
    fn agrees(&self, other: &Self) -> bool {
        self.density.to_bits() == other.density.to_bits()
            && self.sound_speed.to_bits() == other.sound_speed.to_bits()
            && match (self.gas, other.gas) {
                (Some(a), Some(b)) => a == b,
                (None, None) => true,
                _ => false,
            }
    }
}
// The receiver's numerical medium must agree with the complete state being
// attached to the ensemble; a source inlet is not a substitute for exterior gas.
fn receiver_listener(
    receiver: fs_couple::pcm_wav::baffled::RayleighMedium,
    exterior: Option<&GasState>,
) -> Result<Listener, String> {
    if let Some(gas) = exterior {
        if receiver.density.to_bits() != gas.density.to_bits()
            || receiver.sound_speed.to_bits() != gas.sound_speed.to_bits()
        {
            return Err("receiver density/sound speed disagree with the bound exterior gas".into());
        }
    }
    Ok(Listener::Gas(ListenerMedium {
        density: receiver.density,
        sound_speed: receiver.sound_speed,
        gas: exterior.copied(),
    }))
}
// Preserve finite wrappers for producers whose raw scheduler has no horizon.
type DynamicEnsemble = PressureEnsemble<Box<dyn PressureRenderer>>;
fn boxed(renderer: ScheduledRenderer, mut info: SourceInfo) -> (Box<dyn PressureRenderer>, SourceInfo) {
    info.controls += renderer.pending_controls().len();
    (Box::new(renderer), info)
}
fn load(input: &Input, block: usize) -> Result<(Box<dyn PressureRenderer>, SourceInfo), String> {
    load_with_gas_sources(input, block, None)
}
fn load_with_gas_sources(
    input: &Input,
    block: usize,
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(Box<dyn PressureRenderer>, SourceInfo), String> {
    let file = std::fs::File::open(&input.path).map_err(|e| format!("input {}: {e}", input.path.display()))?;
    let max_bytes = match input.kind {
        Kind::Bow => MAX_BOWED_PERFORMANCE_BYTES, Kind::Reed => MAX_REED_PERFORMANCE_BYTES,
        Kind::Flue => MAX_FLUE_PERFORMANCE_BYTES,
        Kind::Lip => MAX_LIP_PERFORMANCE_BYTES, Kind::Voice => MAX_VOICE_PERFORMANCE_BYTES,
        _ => MAX_PART_BYTES,
    };
    let mut bytes = Vec::new();
    file.take((max_bytes + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("input {}: {e}", input.path.display()))?;
    if bytes.len() > max_bytes { return Err(format!("ensemble source exceeds its {max_bytes}-byte input budget")); }
    // A part's material-card packs are named relative to that part's own file.
    let parent = input.path.parent().unwrap_or(std::path::Path::new("."));
    let load_radiation = |name: &str| {
        let file = std::fs::File::open(parent.join(name))
            .map_err(|e| format!("MIMO radiation artifact {name:?} open failed: {e}"))?;
        let mut bytes = Vec::new();
        file.take((MAX_MIMO_RADIATION_BYTES + 1) as u64).read_to_end(&mut bytes)
            .map_err(|e| format!("MIMO radiation artifact {name:?} read failed: {e}"))?;
        if bytes.len() > MAX_MIMO_RADIATION_BYTES {
            return Err("MIMO radiation artifact exceeds its input budget".into());
        }
        Ok(bytes)
    };
    let result = match input.kind {
        Kind::Lip => {
            let description = LipDescription::from_bytes_with_gas_sources(&bytes, sources)?;
            let p = description.build()?;
            let i = p.info();
            let listener = Listener::Gas(ListenerMedium::of(description.ambient_gas()));
            let exterior_json = format!("{},\"observation_scope\":\"lip-driven bore mouth flow through its original baffled exterior receiver\",\"control_item_kind\":\"source pressure/tension initial values and events\",\"mode_count_kind\":\"retained transverse bore order plus one lip coordinate\"",
                super::thermochemical::provenance(std::slice::from_ref(description.gas_record())));
            let info = SourceInfo {
                kind: input.kind, schema: LIP_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.sample_rate_hz, mechanics_samples: i.samples,
                source_filter_delay_samples: p.filter_delay_samples(),
                components: 1, modes: description.transverse_mode_count() + 1,
                controls: description.control_event_count(), exterior_json, listener,
            };
            (Box::new(p) as Box<dyn PressureRenderer>, info)
        }
        Kind::Voice => {
            let description = VoiceDescription::from_bytes_with_gas_sources(&bytes, sources)?;
            let p = description.build()?;
            let i = p.info();
            let listener = Listener::Gas(ListenerMedium::of(description.ambient_gas()));
            let exterior_json = format!("{},\"observation_scope\":\"glottal-driven tract mouth flow through its original multirate baffled exterior receiver\",\"control_item_kind\":\"source subglottal pressure breakpoints\",\"mode_count_kind\":\"one or two glottal fold coordinates\"",
                super::thermochemical::provenance(std::slice::from_ref(description.gas_record())));
            let mechanics_samples = i.samples.checked_mul(u64::from(description.mechanics_per_output))
                .ok_or("voice mechanical window overflow")?;
            let info = SourceInfo {
                kind: input.kind, schema: VOICE_PERFORMANCE_SCHEMA,
                hash: fs_blake3::hash_domain("org.frankensim.fs-couple.voice-performance.v1", &bytes),
                rate: i.output_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.mechanics_rate_hz, mechanics_samples,
                source_filter_delay_samples: p.filter_delay_samples(),
                components: 1, modes: if i.two_mass { 2 } else { 1 },
                controls: description.subglottal.len(), exterior_json, listener,
            };
            (Box::new(p) as Box<dyn PressureRenderer>, info)
        }
        Kind::Modal => {
            let p = ModalPerformance::from_bytes_with_gas_sources(&bytes, block, sources).map_err(|e| e.to_string())?;
            let i = p.info();
            let gas_json = match p.gas_record() {
                Some(record) => super::thermochemical::provenance(std::slice::from_ref(record)),
                None => ",\"radiation\":\"authored, gas-invariant\"".to_string(),
            };
            let listener = p.gas_record().map_or(Listener::Unconstrained("authored transfer; no radiation block"),
                |record| Listener::Gas(ListenerMedium::of(record.state())));
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: i.schema, hash: i.input_hash, rate: i.sample_rate_hz,
                samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.sample_rate_hz, mechanics_samples: i.samples, source_filter_delay_samples: 0.0,
                components: i.voices, modes: i.modes, controls: 0, exterior_json: gas_json, listener,
            })
        }
        Kind::Bow => {
            let p = BowedPerformance::from_bytes_with_gas_sources(&bytes, block, sources).map_err(|e| e.to_string())?;
            let i = p.info();
            // The bridge body carries the gas's baffled piston of its face area.
            let gas_json=format!("{}{}", super::thermochemical::provenance(std::slice::from_ref(p.gas_record())),
                ",\"radiation\":\"baffled piston of the body face area, radius sqrt(A/pi), under the gas record; heard at r/c\"");
            let listener = Listener::Gas(ListenerMedium::of(p.gas_record().state()));
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: BOWED_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.sample_rate_hz, mechanics_samples: i.samples, source_filter_delay_samples: 0.0,
                components: 1, modes: i.string_modes + 1, controls: i.compiled_controls, exterior_json: gas_json, listener,
            })
        }
        Kind::Plate => {
            let p = PlatePerformance::from_bytes_with_all_sources(&bytes, block, sources,
                super::plate_input::section_pack_loader(parent)).map_err(|e| e.to_string())?;
            let i = p.info();
            let gas_json=format!("{}{}", super::thermochemical::provenance(std::slice::from_ref(p.gas_record())),
                super::plate_input::section_sources_json(p.section_provenance()));
            let listener = Listener::Gas(ListenerMedium::of(p.gas_record().state()));
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: PLATE_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.sample_rate_hz, mechanics_samples: i.samples, source_filter_delay_samples: 0.0,
                components: 1, modes: i.retained_modes, controls: 0, exterior_json: gas_json, listener,
            })
        }
        Kind::Flue => match std::str::from_utf8(&bytes)
            .map_err(|_| "flue performance must be UTF-8")?
            .lines().next().unwrap_or("").trim()
        {
            FLUE_PERFORMANCE_SCHEMA => {
                let p = FlueFilePerformance::from_bytes_with_all_sources(
                    &bytes, block, &CancelGate::new(), sources,
                    super::plate_input::section_pack_loader(parent), load_radiation,
                ).map_err(|e| e.to_string())?;
                let i = p.info();
                let renderer = p.renderer();
                // The source owns the complete phrase, all mechanical substeps
                // and its differentiating microphone. Mix its output pressure.
                let exterior = renderer.radiation_input()
                    .ok_or("flue microphone requires its declared complete exterior gas")?;
                let listener = Listener::Gas(ListenerMedium::of(&exterior.exterior));
                let source_filter_delay_samples = renderer.receiver().filter_delay_samples();
                let mechanics_samples = i.samples.checked_mul(renderer.mechanics_per_output() as u64)
                    .ok_or("flue mechanical window overflow")?;
                let exterior_json = format!(",{},\"observation_scope\":\"exterior pressure from the original multirate flue receiver\"",
                    super::wind_input::flue_provenance(&p));
                (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                    kind: input.kind, schema: FLUE_PERFORMANCE_SCHEMA, hash: i.input_hash,
                    rate: i.sample_rate_hz, samples: i.samples,
                    mechanics_rate: i.mechanics_rate_hz, mechanics_samples,
                    source_filter_delay_samples, source_full_scale_pa: i.full_scale_pa,
                    components: 1, modes: 0, controls: i.compiled_controls, exterior_json, listener,
                })
            }
            FLUE_DG_PERFORMANCE_SCHEMA => {
                let p = FlueDgFilePerformance::from_bytes_with_sources(
                    &bytes, block, &CancelGate::new(), sources,
                    super::plate_input::section_pack_loader(parent),
                ).map_err(|e| e.to_string())?;
                let i = p.info();
                let renderer = p.renderer();
                let microphone = renderer.baffled_receiver();
                let listener = if microphone.is_some() {
                    let exterior = renderer.voice().dg()
                        .ok_or("DG flue microphone lacks its original bore")?
                        .load().radiation_record()
                        .ok_or("DG flue microphone lacks its original computed radiation source")?
                        .gas_record().state();
                    Listener::Gas(ListenerMedium::of(exterior))
                } else {
                    Listener::Unconstrained("internal DG inlet trace; no exterior microphone")
                };
                let source_filter_delay_samples = microphone.map_or(0.0, |m| m.filter_delay_samples());
                let scope = if microphone.is_some() {
                    "computed DG flue terminal and optional mouth flows through their original baffled exterior receivers"
                } else { "internal DG step-mean inlet pressure; not an exterior microphone" };
                let exterior_json = format!(",{},\"observation_scope\":{}",
                    super::wind_input::flue_dg_provenance(&p), json_string(scope));
                (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                    kind: input.kind, schema: FLUE_DG_PERFORMANCE_SCHEMA, hash: i.input_hash,
                    rate: i.sample_rate_hz, samples: i.samples,
                    mechanics_rate: i.mechanics_rate_hz, mechanics_samples: i.samples,
                    source_filter_delay_samples, source_full_scale_pa: i.full_scale_pa,
                    components: 1, modes: 0, controls: i.compiled_controls, exterior_json, listener,
                })
            }
            _ => return Err("--flue requires a supported flue or flue-DG source schema".into()),
        },
        Kind::Valve => match std::str::from_utf8(&bytes)
            .map_err(|_| "valve performance must be UTF-8")?
            .lines().next().unwrap_or("").trim()
        {
            MATERIAL_MODAL_DG_VALVE_SCHEMA => {
                let p = MaterialModalDgValvePerformance::from_bytes_with_sources(
                    &bytes, block, &CancelGate::new(), sources,
                    super::plate_input::section_pack_loader(parent),
                ).map_err(|e| e.to_string())?;
                let i = p.info();
                let renderer = p.renderer();
                let microphone = renderer.baffled_receiver();
                let listener = if microphone.is_some() {
                    let exterior = renderer.system().load().radiation_record()
                        .ok_or("DG microphone lacks its original computed radiation source")?
                        .gas_record().state();
                    Listener::Gas(ListenerMedium::of(exterior))
                } else {
                    Listener::Unconstrained("internal DG inlet trace; no exterior microphone")
                };
                let source_filter_delay_samples = microphone.map_or(0.0, |m| m.filter_delay_samples());
                let controls = renderer.pending_controls().len() + (0..i.blades).map(|blade|
                    renderer.pending_plate_force_controls(blade).len()
                        + renderer.pending_material_pressure_controls(blade).len()).sum::<usize>();
                let scope = if microphone.is_some() {
                    "computed DG terminal flow through its original baffled exterior receiver"
                } else { "internal DG step-mean inlet pressure; not an exterior microphone" };
                let exterior_json = format!(",{},\"observation_scope\":{}",
                    super::wind_input::material_modal_dg_provenance(&p), json_string(scope));
                (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                    kind: input.kind, schema: MATERIAL_MODAL_DG_VALVE_SCHEMA, hash: i.input_hash,
                    rate: i.output_rate_hz, samples: i.samples,
                    mechanics_rate: i.mechanics_rate_hz, mechanics_samples: i.samples,
                    source_filter_delay_samples, source_full_scale_pa: i.full_scale_pa,
                    components: i.blades, modes: i.modes, controls, exterior_json, listener,
                })
            }
            MATERIAL_MODAL_VALVE_SCHEMA => {
                let p = MaterialModalValvePerformance::from_bytes_with_all_sources(
                    &bytes, block, &CancelGate::new(), sources,
                    super::plate_input::section_pack_loader(parent), load_radiation,
                ).map_err(|e| e.to_string())?;
                let i = p.info();
                let renderer = p.renderer();
                // This source already owns every mechanical substep and its
                // differentiating microphone. Ensemble receives OUTPUT frames.
                let exterior = renderer.source().radiation().exterior.as_ref()
                    .ok_or("material-modal microphone requires its declared complete exterior gas")?;
                let listener = Listener::Gas(ListenerMedium::of(exterior));
                let source_filter_delay_samples = renderer.source().receiver().filter_delay_samples();
                let mechanics_samples = i.samples.checked_mul(renderer.source().mechanics_per_output() as u64)
                    .ok_or("material-modal mechanical window overflow")?;
                let controls = renderer.pending_controls().len()
                    + [ValveSide::Upstream, ValveSide::Downstream].into_iter().map(|side|
                        renderer.pending_boundary_controls(side).len()).sum::<usize>()
                    + (0..i.blades).map(|blade| renderer.pending_plate_force_controls(blade).len()
                        + renderer.pending_material_pressure_controls(blade).len()).sum::<usize>();
                let exterior_json = format!(",{},\"observation_scope\":\"exterior pressure from the original multirate material-modal receiver\"",
                    super::wind_input::material_modal_provenance(&p));
                (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                    kind: input.kind, schema: MATERIAL_MODAL_VALVE_SCHEMA, hash: i.input_hash,
                    rate: i.output_rate_hz, samples: i.samples,
                    mechanics_rate: i.mechanics_rate_hz, mechanics_samples, source_filter_delay_samples,
                    source_full_scale_pa: i.full_scale_pa, components: i.blades, modes: i.modes,
                    controls, exterior_json, listener,
                })
            }
            PLATE_VALVE_PERFORMANCE_SCHEMA => {
            let p = PlateValvePerformance::from_bytes_with_all_sources(&bytes, block, &CancelGate::new(),sources,
                super::plate_input::section_pack_loader(parent)).map_err(|e|e.to_string())?;
            let i = p.info();
            let exterior_json = format!("{}{}{}{}",super::wind_input::outlet_provenance(p.renderer(),i.radiation_load),
                super::duct_metadata::provenance(&p),super::thermochemical::provenance(p.gas_records()),
                super::plate_input::section_sources_json(p.section_provenance()));
            // An active microphone uses the same actual bound exterior as its
            // reaction/observer; without one it hears its section's sourced gas.
            // Internal traces retain the explicit-file policy.
            let listener = match p.renderer().baffled_medium() {
                Some(medium) => receiver_listener(
                    medium, super::duct_metadata::aperture_exterior_gas(&p)
                        .or_else(|| p.renderer().baffled_gas()),
                )?,
                None => p.exterior_gas().map_or(
                    Listener::Unconstrained("internal trace; no declared exterior"),
                    |exterior| Listener::Gas(ListenerMedium::of(exterior)),
                ),
            };
            let source_filter_delay_samples = p.renderer().baffled_receiver()
                .map_or(0.0, |receiver| receiver.filter_delay_samples());
            (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                kind: input.kind, schema: PLATE_VALVE_PERFORMANCE_SCHEMA, hash:i.input_hash,
                rate:i.sample_rate_hz,samples:i.samples,source_full_scale_pa:i.full_scale_pa,
                mechanics_rate:i.sample_rate_hz,mechanics_samples:i.samples,source_filter_delay_samples,
                components:1,modes:1,controls:i.compiled_controls+i.force_controls,exterior_json,listener,
            })
            }
            _ => return Err("--valve requires a supported plate-valve, material-modal or material-modal-DG source schema".into()),
        },
        Kind::Reed => {
            let p = ReedPerformance::from_bytes_with_gas_sources(&bytes, block,
                fs_couple::pcm_wav::observation::ObservationPolicy::Diagnostic,sources).map_err(|e| e.to_string())?;
            let i = p.info();
            let exterior_json = format!("{}{}",super::duct_metadata::reed_provenance(&p),super::thermochemical::provenance(p.gas_records()));
            // The baffled outlet's one bound exterior; the legacy proxy has none.
            let listener = p.outlet_receiver().map_or(Listener::Unconstrained("legacy listener_m bore/jet proxy"),
                |o| Listener::Gas(ListenerMedium::of(&o.gas)));
            (Box::new(p) as Box<dyn PressureRenderer>, SourceInfo {
                kind: input.kind, schema: REED_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                mechanics_rate: i.sample_rate_hz, mechanics_samples: i.samples, source_filter_delay_samples: 0.0,
                components: 1, modes: 0, controls: i.compiled_controls, exterior_json, listener,
            })
        }
    };
    Ok(result)
}

fn prepare(options: &Options) -> Result<(DynamicEnsemble, Vec<SourceInfo>, u64), String> {
    prepare_with_gas_sources(options, None)
}
fn prepare_with_gas_sources(
    options: &Options,
    gas_sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(DynamicEnsemble, Vec<SourceInfo>, u64), String> {
    let mut parts = Vec::new();
    let mut sources = Vec::new();
    let mut source_delays = Vec::new();
    let mut horizon = None;
    let mut converted = false;
    let mut listener: Option<(usize, Kind, ListenerMedium)> = None;
    let (mut components, mut modes, mut controls) = (0_usize, 0_usize, 0_usize);
    for (index, input) in options.inputs.iter().enumerate() {
        let (renderer, info) = load_with_gas_sources(input, options.block,gas_sources)
            .map_err(|e| format!("ensemble part {index} ({}): {e}", input.kind.label()))?;
        if info.rate != RATE && !options.decimate {
            return Err(format!("ensemble part {index} needs explicit --decimate; its pressure producer outputs {} Hz", info.rate));
        }
        converted |= info.rate != RATE;
        // One observer sits in one gas: every part that declares the gas it is
        // heard in must declare the same one.
        if let Listener::Gas(medium) = info.listener {
            match listener {
                None => listener = Some((index, info.kind, medium)),
                Some((first, kind, shared)) if !shared.agrees(&medium) => return Err(format!(
                    "ensemble part {index} ({}) is heard in a different gas than part {first} ({}): rho {:e} kg/m3, c {:e} m/s against rho {:e} kg/m3, c {:e} m/s; one ensemble has one listener gas",
                    info.kind.label(), kind.label(), medium.density, medium.sound_speed, shared.density, shared.sound_speed)),
                Some(_) => {}
            }
        }
        components = components.checked_add(info.components).ok_or("ensemble component count overflow")?;
        modes = modes.checked_add(info.modes).ok_or("ensemble mode count overflow")?;
        controls = controls.checked_add(info.controls).ok_or("ensemble control count overflow")?;
        if components > MAX_COMPONENTS || modes > MAX_MODES || controls > MAX_CONTROLS {
            return Err("ensemble exceeds 128 source components, 8192 modes or 262144 compiled controls".into());
        }
        let observed = DecimatedRenderer::new(renderer, info.rate, RATE, options.block)
            .map_err(|e| format!("ensemble part {index} clock: {e}"))?;
        // Receiver latency is already present in its pressure trace. Convert
        // only its units; the shared ensemble owner validates exact integrality
        // and adds the remaining pure delay. Physical r/c stays untouched.
        source_delays.push(info.source_filter_delay_samples / observed.info().ratio as f64);
        let samples = observed.output_samples_for(info.samples).map_err(|e| e.to_string())?;
        if samples == 0 || samples > MAX_SAMPLES { return Err("ensemble duration must be in (0,600] seconds".into()); }
        if let Some(expected) = horizon {
            if samples != expected { return Err("all source files must declare exactly the same physical duration; no truncation or padding".into()); }
        } else { horizon = Some(samples); }
        sources.push(info);
        parts.push(observed);
    }
    if options.decimate && !converted {
        return Err("--decimate requires at least one higher-rate source; 48000-Hz parts need no conversion".into());
    }
    let samples = horizon.ok_or_else(|| "ensemble has no source duration".to_string())?;
    let renderer = PressureEnsemble::new_with_source_delays(parts, PressureEnsembleConfig {
        sample_rate_hz: RATE, samples, max_block: options.block, max_parts: MAX_PARTS,
    }, &source_delays).map_err(|e| e.to_string())?;
    renderer.validate_sample_count(samples).map_err(|e| e.to_string())?;
    Ok((renderer, sources, samples))
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    run_with_gas_sources(args, None)
}
pub(super) fn run_with_gas_sources(
    args: &[String],
    gas_sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    let options = options(args)?;
    let output = &options.output;
    let sidecar = output.with_extension("provenance.json");
    if output == &sidecar || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let (mut renderer, sources, samples) = prepare_with_gas_sources(&options, gas_sources)?;
    let sample_count = usize::try_from(samples).map_err(|_| "output sample count exceeds this host".to_string())?;
    // Finish every source admission and reduction before creating either file.
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, sample_count, options.block, options.full_scale_pa,
    )?;
    let parts: Vec<_> = sources.iter().zip(renderer.part_info()).map(|(source, aligned)| {
        let o = aligned.observation;
        let mechanics_rate = source.mechanics_rate;
        let mechanics_samples = source.mechanics_samples;
        let source_delay = aligned.source_filter_delay_output_samples;
        let total_delay = aligned.total_filter_delay_output_samples;
        let default_scope = if source.exterior_json.contains("\"observation_scope\":") {
            ""
        } else if source.kind == Kind::Reed {
            ",\"observation_scope\":\"bore-pressure-plus-compact-jet-proxy; not exterior microphone\""
        } else if source.kind == Kind::Valve {
            ",\"observation_scope\":\"internal coupled tube pressure; not an exterior microphone\""
        } else { "" };
        let scope = source.exterior_json.as_str();
        format!("{{\"kind\":\"{}\",\"schema\":\"{}\",\"blake3\":\"{}\",\
            \"mechanics_sample_rate_hz\":{mechanics_rate},\"mechanics_samples\":{mechanics_samples},\
            \"pressure_sample_rate_hz\":{},\"pressure_samples\":{},\"ratio\":{},\
            \"filter\":\"{}\",\"filter_delay_output_samples\":{},\"alignment_delay_samples\":{},\
            \"source_filter_delay_output_samples\":{source_delay},\"total_filter_delay_output_samples\":{total_delay},\
            \"first_output_source_index\":{},\"source_full_scale_pa\":{:e},\"source_pcm_scale_applied\":false,\
            \"components\":{},\"modes\":{},\"compiled_controls\":{}{default_scope}{scope}}}",
            source.kind.label(), source.schema, source.hash.to_hex(), source.rate, source.samples,
            o.ratio, o.filter_profile, o.delay_output_samples, aligned.alignment_delay_samples,
            o.first_output_source_index, source.source_full_scale_pa, source.components, source.modes, source.controls)
    }).collect();
    let hash = rendered.hash.to_hex();
    let unconstrained: Vec<_> = sources.iter().enumerate().filter_map(|(i, s)| match s.listener {
        Listener::Unconstrained(reason) => Some(format!("{{\"part\":{i},\"kind\":\"{}\",\"reason\":\"{reason}\"}}", s.kind.label())),
        Listener::Gas(_) => None,
    }).collect();
    let listener = match sources.iter().enumerate().find_map(|(i, s)| match s.listener {
        Listener::Gas(m) => Some((i, m)), Listener::Unconstrained(_) => None,
    }) {
        Some((first, m)) => format!("{{\"medium\":{},\"set_by_part\":{first},\"unconstrained_parts\":[{}]}}",
            m.provenance_json(), unconstrained.join(",")),
        None => format!("{{\"declared\":\"none\",\"unconstrained_parts\":[{}]}}", unconstrained.join(",")),
    };
    let provenance = format!("{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"ensemble-input\",\
        \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{},\"full_scale_pa\":{:e},\
        \"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{hash}\",\
        \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
        \"ensemble\":{{\"sum_order\":\"input-order\",\"common_delay_output_samples\":{},\
        \"delay_compensated\":false,\"alignment_scope\":\"declared receiver numerical filter plus outer decimator; physical propagation retained\",\
        \"initial_history\":\"zero\",\"tail\":\"no-flush-declared-window\",\
        \"listener_medium\":{listener},\"parts\":[{}],\"model_scope\":\"independent source pressures; common observer and time origin required; no new mechanical coupling or physical-validation claim\"}}}}",
        options.block, options.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        renderer.delay_output_samples(), parts.join(","));
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!("{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"ensemble-input\",\
        \"wav\":{},\"samples\":{samples},\"parts\":{},\"clipped\":{},\"peak_pa\":{:e},\
        \"rms_pa\":{:e},\"wav_blake3\":\"{hash}\"}}",
        json_string(&output.display().to_string()), sources.len(), rendered.clipped, rendered.peak_pa, rendered.rms_pa);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_string()).collect() }
    #[test]
    fn observed_medium_cannot_be_receipted_as_a_different_source_gas() {
        let inlet = GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap();
        let exterior = GasState::try_new_moist_air(293.15, 81060.0, 0.0).unwrap();
        let observed = fs_couple::pcm_wav::baffled::RayleighMedium::from(&exterior);
        assert_eq!(
            receiver_listener(observed, Some(&exterior)).unwrap(),
            Listener::Gas(ListenerMedium::of(&exterior)),
        );
        assert!(receiver_listener(observed, Some(&inlet)).is_err());
        assert!(!ListenerMedium::of(&exterior).agrees(&ListenerMedium::of(&inlet)));
        let receiver_only = match receiver_listener(observed, None).unwrap() {
            Listener::Gas(medium) => medium,
            Listener::Unconstrained(_) => unreachable!("receiver declares a gas medium"),
        };
        assert_eq!(
            receiver_only,
            ListenerMedium {
                density: exterior.density, sound_speed: exterior.sound_speed, gas: None,
            },
        );
        assert!(
            !receiver_only.agrees(&ListenerMedium::of(&exterior)),
            "matching rho/c without a retained complete gas cannot prove material compatibility"
        );
        assert!(
            receiver_only.agrees(&receiver_only),
            "two receiver-only media with the same consumed coefficients remain compatible"
        );
    }
    #[test]
    fn listener_provenance_distinguishes_complete_state_from_receiver_coefficients() {
        let gas = GasState::try_new_moist_air(293.15, 101325.0, 0.35).unwrap();
        let complete = ListenerMedium::of(&gas).provenance_json();
        assert!(complete.contains("\"state_authority\":\"complete_gas_state\""));
        assert!(complete.contains("\"temperature_k\":"));
        assert!(complete.contains("\"pressure_pa\":"));
        assert!(complete.contains("\"water_mole_fraction\":"));

        let receiver = ListenerMedium {
            density: gas.density,
            sound_speed: gas.sound_speed,
            gas: None,
        }
        .provenance_json();
        assert!(receiver.contains("\"state_authority\":\"receiver_coefficients_only\""));
        assert!(!receiver.contains("\"temperature_k\":"));
        assert!(!receiver.contains("\"pressure_pa\":"));
    }

    #[test]
    fn physical_duration_stays_in_sources_and_pcm_scale_is_explicit() {
        let valid = args(&["out.wav", "--plate", "a.performance", "--modal", "b.performance", "--full-scale-pa", "2", "--block", "37", "--decimate"]);
        let o = options(&valid).unwrap();
        assert_eq!(o.inputs.iter().map(|i| i.kind).collect::<Vec<_>>(), vec![Kind::Plate, Kind::Modal]);
        assert_eq!(o.full_scale_pa, 2.0); assert_eq!(o.block, 37); assert!(o.decimate);
        for flag in ["--seconds", "--samples", "--schedule", "--gain", "--unknown", "--decimate"] {
            let mut a = valid.clone(); a.push(flag.into());
            assert!(options(&a).is_err());
        }
        assert!(options(&args(&["out.wav", "--modal", "a.performance"])).is_err());
        assert!(options(&args(&["out.wav", "--full-scale-pa", "2"])).is_err());
        for value in ["0", "-1", "NaN", "inf"] {
            assert!(options(&args(&["out.wav", "--modal", "a", "--full-scale-pa", value])).is_err());
        }
    }
}
