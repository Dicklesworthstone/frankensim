//! Supplied primitive or geometry-derived valve pressure to the sole WAV owner.
use std::io::{Read, Write};
use std::path::Path;
use fs_blake3::ContentHash;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::render::schedule::reed::{ReedPerformance, REED_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::performance::{ApertureObservation, AperturePerformance, CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::{PlateValvePerformance, MAX_PLATE_VALVE_PERFORMANCE_BYTES, PLATE_VALVE_PERFORMANCE_SCHEMA};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_exec::CancelGate;
use super::{RATE, create_outputs, json_string, stream_output};

fn options(args: &[String]) -> Result<(&str, &str, usize, bool), String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut decimate = false;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--decimate" => {
                if decimate { return Err("--decimate may only be supplied once".into()); }
                decimate = true;
            }
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                block = Some(iter.next().and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?);
            }
            value if value.starts_with('-') => return Err(format!(
                "unsupported wind option {value:?}; geometry, pressure, duration and full-scale belong to the source file"
            )),
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render wind INPUT.performance OUT.wav [--block N] [--decimate]".into());
    };
    Ok((*input, *output, block.unwrap_or(512), decimate))
}

// A physical propagation delay is source geometry, not an output/filter latency
// to be removed by ensemble alignment. Keep it on the mechanical sample clock.
pub(super) fn outlet_provenance(p: &AperturePerformance, radiation: Option<BaffledRadiationLoad>) -> String {
    let load_json=if let Some(load)=radiation {
        let spec=load.load();let term=spec.terms()[0];
        format!(",\"radiation_load\":{{\"model\":\"compact-baffled-piston-positive-real-v1\",\"maximum_frequency_hz\":{:e},\"resistance_pa_s_m3\":{:e},\"pole_rate_per_s\":{:e},\"checked_max_complex_relative_error\":{:e},\"checked_max_resistance_relative_error\":{:e},\"replaces_memoryless_reflection\":true,\"extra_end_correction_m\":0,\"scope\":\"finite-sample low-ka analogue and bilinear checks, not broadband or measured calibration\"}}",
            load.maximum_frequency_hz(),term.resistance_pa_s_m3,term.rate_per_s,
            load.max_complex_relative_error(),load.max_resistance_relative_error())
    } else {String::new()};
    let (c,radius,speed)=match (p.observation(),p.system()) {
        (ApertureObservation::TubeBaffled(c),CoupledAperture::Tube(t))=>(c,t.spec().radius_m,t.spec().sound_speed_m_s),
        (ApertureObservation::NetworkBaffled {node,receiver},CoupledAperture::Network(n))=> {
            let section=n.spec().sections.iter().find(|s|s.nodes.contains(&node)).expect("admitted physical outlet section");
            (receiver,section.radius_m,n.spec().sound_speed_m_s)
        }
        _=>return if radiation.is_some() {
            format!(",\"observation_scope\":\"internal coupled tube pressure; not an exterior microphone\"{load_json}")
        } else {String::new()},
    };
    let mic=p.baffled_receiver().expect("admitted receiver");
    let scope=if radiation.is_some() {
        "exterior Rayleigh pressure from radiation-loaded terminal flow; compact passive load and uniform circular outlet in infinite rigid baffle"
    } else {
        "one-way exterior Rayleigh pressure; uniform circular outlet in infinite rigid baffle; terminal load is independently supplied"
    };
    format!(",\"observation_scope\":\"{scope}\",\"outlet_receiver\":{{\"position_m\":[{:e},{:e},{:e}],\"radius_m\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"propagation_delay_mechanical_samples\":[{},{}],\"radiation_feedback_added\":{}}}{load_json}",
        c.position_m[0],c.position_m[1],c.position_m[2],radius,
        p.system().aperture().spec().density_kg_m3,speed,
        c.radial_rings,c.angular_points,c.maximum_frequency_hz,mic.delay_samples.0,mic.delay_samples.1,radiation.is_some())
}

struct Loaded {
    source: Box<dyn PressureRenderer>, rate:u32, samples:u64, full_scale_pa:f64,
    hash:ContentHash, fixture:&'static str, source_json:String,
}
fn load(bytes:&[u8],block:usize)->Result<Loaded,String> {
    // Exact versioned source header selects the existing owner; filenames and
    // instrument labels never choose geometry or material constants.
    if bytes.starts_with(b"frankensim-plate-valve-performance-v1\n") {
        let p=PlateValvePerformance::from_bytes(bytes,block,&CancelGate::new()).map_err(|e|e.to_string())?;
        let i=p.info();let a=p.renderer().system().aperture();
        let plate=a.plate_reduction().expect("plate file retains its specimen");
        let observation=match p.renderer().observation() {
            ApertureObservation::Inlet=>"inlet",
            ApertureObservation::TubeTerminal=>"terminal",
            ApertureObservation::NetworkNode(node)=>match p.renderer().system() {
                CoupledAperture::Network(n) if i.duct_sections==1 && matches!(n.spec().nodes[node],
                    fs_couple::bernoulli_aperture::network::NetworkNode::Termination {..}
                    | fs_couple::bernoulli_aperture::network::NetworkNode::Impedance {..}
                    | fs_couple::bernoulli_aperture::network::NetworkNode::Relaxation {..})=>"terminal",
                _=>"network-node",
            },
            ApertureObservation::TubeBaffled(_)|ApertureObservation::NetworkBaffled {..}=>"baffled-outlet",
        };
        let receiver_json=format!("{}{}",outlet_provenance(p.renderer(),i.radiation_load),super::duct_metadata::provenance(&p));
        let scope=if observation!="baffled-outlet" { "internal pressure, not an exterior microphone" }
            else if i.radiation_load.is_some() { "exterior baffled-outlet pressure with compact passive radiation feedback; not broadband matched radiation or measured calibration" }
            else { "one-way exterior baffled-outlet pressure; no matched radiation load or measured-instrument claim" };
        let requested=match p.renderer().system() {
            CoupledAperture::Tube(t)=>t.spec().length_m,
            CoupledAperture::Network(n)=>n.spec().sections.iter().map(|s|s.length_m).sum(),
        };
        let source_json=format!("\"plate_valve_input\":{{\"schema\":\"{PLATE_VALVE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{}\",\"nodes\":{},\"triangles\":{},\"sections\":{},\"memory_branches\":{},\"compiled_controls\":{},\"pressure_point\":\"{observation}\",\"requested_tube_length_m\":{:e},\"represented_tube_length_m\":{:e},\"effective_mass_kg\":{:e},\"stiffness_n_m\":{:e},\"pressure_area_m2\":{:e},\"model_scope\":\"one linear plate mode; spatial lay and supplied material history; lossless propagating sections with explicit local loads; {scope}\"{receiver_json}}}",
            i.input_hash.to_hex(),i.nodes,i.triangles,i.sections,i.memory_branches,i.compiled_controls,
            requested,i.represented_tube_length_m,plate.mass_kg(),plate.stiffness_n_m(),plate.pressure_area_m2());
        Ok(Loaded{source:Box::new(p.into_renderer()),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"plate-valve-input",source_json})
    } else {
        let p=ReedPerformance::from_bytes(bytes,block).map_err(|e|e.to_string())?;
        let i=p.info();
        let source_json=format!("\"reed_input\":{{\"schema\":\"{REED_PERFORMANCE_SCHEMA}\",\"blake3\":\"{}\",\"segments\":{},\"tone_holes\":{},\"massive_reed\":{},\"gesture_events\":{},\"compiled_controls\":{},\"model_scope\":\"authored reed primitives and static duct; bore-pressure plus compact-jet proxy; not a calibrated exterior microphone or material-identification claim\"}}",
            i.input_hash.to_hex(),i.segments,i.tone_holes,i.massive_reed,i.gesture_events,i.compiled_controls);
        Ok(Loaded{source:Box::new(p),rate:i.sample_rate_hz,samples:i.samples,
            full_scale_pa:i.full_scale_pa,hash:i.input_hash,fixture:"reed-input",source_json})
    }
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    let (input, output, block, decimate) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("wind input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_PLATE_VALVE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("wind input read failed: {e}"))?;
    let loaded=load(&bytes,block)?;
    if !decimate && loaded.rate != RATE {
        return Err("non-48000-Hz mechanics require explicit --decimate; no implicit resampling".into());
    }
    if decimate && loaded.rate <= RATE {
        return Err("--decimate requires mechanics above 48000 Hz; ordinary 48000-Hz input needs no conversion".into());
    }
    // The entire finite physical runtime stays inside this observer, including
    // spatial contact and material memory. Only its output pressure is filtered.
    let mut renderer = DecimatedRenderer::new(loaded.source, loaded.rate, RATE, block)
        .map_err(|e| e.to_string())?;
    let samples = renderer.output_samples_for(loaded.samples).map_err(|e| e.to_string())?;
    if samples == 0 { return Err("wind source needs at least one complete output interval".into()); }
    renderer.validate_sample_count(samples).map_err(|e| e.to_string())?;
    let samples = usize::try_from(samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let observation = renderer.info();
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, samples, block, loaded.full_scale_pa,
    )?;
    let input_hash = loaded.hash.to_hex();let wav_hash = rendered.hash.to_hex();let fixture=loaded.fixture;
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"{fixture}\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",{},\
         \"observation\":{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\
         \"output_sample_rate_hz\":{RATE},\"ratio\":{},\"filter\":\"{}\",\
         \"delay_output_samples\":{},\"first_output_source_index\":{},\
         \"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}}}",
        loaded.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,loaded.source_json,
        loaded.rate, loaded.samples, observation.ratio, observation.filter_profile,
        observation.delay_output_samples, observation.first_output_source_index,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!("{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"{fixture}\",\
        \"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
        \"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa, rendered.rms_pa);
    Ok(())
}
