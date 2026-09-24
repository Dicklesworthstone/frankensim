//! File-driven multi-rate physical ensembles; use existing model and PCM owners.
use std::io::{Read, Write};
use std::path::PathBuf;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_couple::bernoulli_aperture::performance::file::{PlateValvePerformance, PLATE_VALVE_PERFORMANCE_SCHEMA};
use fs_couple::bowed_string::runtime::schedule::file::{BowedPerformance, BOWED_PERFORMANCE_SCHEMA, MAX_BOWED_PERFORMANCE_BYTES};
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
enum Kind { Modal, Plate, Bow, Reed, Valve }
impl Kind {
    fn label(self) -> &'static str {
        match self { Self::Modal => "modal", Self::Plate => "plate", Self::Bow => "bow", Self::Reed => "reed", Self::Valve => "valve" }
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
            "--modal" | "--plate" | "--bow" | "--reed" | "--valve" => {
                if inputs.len() == MAX_PARTS { return Err("ensemble exceeds the 16-part budget".into()); }
                let path = iter.next().filter(|s| !s.starts_with("--"))
                    .ok_or_else(|| format!("{arg} requires a source performance path"))?;
                inputs.push(Input { kind: match arg.as_str() { "--modal" => Kind::Modal, "--plate" => Kind::Plate, "--bow" => Kind::Bow, "--reed" => Kind::Reed, _ => Kind::Valve },
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
                if output.is_some() { return Err("ensemble accepts exactly one output path; inputs require --modal/--plate/--bow/--reed/--valve".into()); }
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
    rate: u32,
    samples: u64,
    source_full_scale_pa: f64,
    components: usize,
    modes: usize,
    controls: usize,
    exterior_json: String,
}
// Preserve finite wrappers for producers whose raw scheduler has no horizon.
type DynamicEnsemble = PressureEnsemble<Box<dyn PressureRenderer>>;
fn boxed(renderer: ScheduledRenderer, mut info: SourceInfo) -> (Box<dyn PressureRenderer>, SourceInfo) {
    info.controls += renderer.pending_controls().len();
    (Box::new(renderer), info)
}
fn load(input: &Input, block: usize) -> Result<(Box<dyn PressureRenderer>, SourceInfo), String> {
    let file = std::fs::File::open(&input.path).map_err(|e| format!("input {}: {e}", input.path.display()))?;
    let max_bytes = match input.kind {
        Kind::Bow => MAX_BOWED_PERFORMANCE_BYTES, Kind::Reed => MAX_REED_PERFORMANCE_BYTES,
        _ => MAX_PART_BYTES,
    };
    let mut bytes = Vec::new();
    file.take((max_bytes + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("input {}: {e}", input.path.display()))?;
    if bytes.len() > max_bytes { return Err(format!("ensemble source exceeds its {max_bytes}-byte input budget")); }
    let result = match input.kind {
        Kind::Modal => {
            let p = ModalPerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
            let i = p.info();
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: i.schema, hash: i.input_hash, rate: i.sample_rate_hz,
                samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                components: i.voices, modes: i.modes, controls: 0, exterior_json: String::new(),
            })
        }
        Kind::Bow => {
            let p = BowedPerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
            let i = p.info();
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: BOWED_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                components: 1, modes: i.string_modes + 1, controls: i.compiled_controls, exterior_json: String::new(),
            })
        }
        Kind::Plate => {
            let p = PlatePerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
            let i = p.info();
            boxed(p.into_renderer(), SourceInfo {
                kind: input.kind, schema: PLATE_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                components: 1, modes: i.retained_modes, controls: 0, exterior_json: String::new(),
            })
        }
        Kind::Valve => {
            let p = PlateValvePerformance::from_bytes(&bytes, block, &CancelGate::new()).map_err(|e|e.to_string())?;
            let i = p.info();
            let exterior_json = super::wind_input::outlet_provenance(p.renderer());
            (Box::new(p.into_renderer()) as Box<dyn PressureRenderer>, SourceInfo {
                kind: input.kind, schema: PLATE_VALVE_PERFORMANCE_SCHEMA, hash:i.input_hash,
                rate:i.sample_rate_hz,samples:i.samples,source_full_scale_pa:i.full_scale_pa,
                components:1,modes:1,controls:i.compiled_controls,exterior_json,
            })
        }
        Kind::Reed => {
            let p = ReedPerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
            let i = p.info();
            (Box::new(p) as Box<dyn PressureRenderer>, SourceInfo {
                kind: input.kind, schema: REED_PERFORMANCE_SCHEMA, hash: i.input_hash,
                rate: i.sample_rate_hz, samples: i.samples, source_full_scale_pa: i.full_scale_pa,
                components: 1, modes: 0, controls: i.compiled_controls, exterior_json: String::new(),
            })
        }
    };
    Ok(result)
}

fn prepare(options: &Options) -> Result<(DynamicEnsemble, Vec<SourceInfo>, u64), String> {
    let mut parts = Vec::new();
    let mut sources = Vec::new();
    let mut horizon = None;
    let mut converted = false;
    let (mut components, mut modes, mut controls) = (0_usize, 0_usize, 0_usize);
    for (index, input) in options.inputs.iter().enumerate() {
        let (renderer, info) = load(input, options.block)
            .map_err(|e| format!("ensemble part {index} ({}): {e}", input.kind.label()))?;
        if info.rate != RATE && !options.decimate {
            return Err(format!("ensemble part {index} needs explicit --decimate; mechanics run at {} Hz", info.rate));
        }
        converted |= info.rate != RATE;
        components = components.checked_add(info.components).ok_or("ensemble component count overflow")?;
        modes = modes.checked_add(info.modes).ok_or("ensemble mode count overflow")?;
        controls = controls.checked_add(info.controls).ok_or("ensemble control count overflow")?;
        if components > MAX_COMPONENTS || modes > MAX_MODES || controls > MAX_CONTROLS {
            return Err("ensemble exceeds 128 source components, 8192 modes or 262144 compiled controls".into());
        }
        let observed = DecimatedRenderer::new(renderer, info.rate, RATE, options.block)
            .map_err(|e| format!("ensemble part {index} clock: {e}"))?;
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
    let renderer = PressureEnsemble::new(parts, PressureEnsembleConfig {
        sample_rate_hz: RATE, samples, max_block: options.block, max_parts: MAX_PARTS,
    }).map_err(|e| e.to_string())?;
    renderer.validate_sample_count(samples).map_err(|e| e.to_string())?;
    Ok((renderer, sources, samples))
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    let options = options(args)?;
    let output = &options.output;
    let sidecar = output.with_extension("provenance.json");
    if output == &sidecar || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let (mut renderer, sources, samples) = prepare(&options)?;
    let sample_count = usize::try_from(samples).map_err(|_| "output sample count exceeds this host".to_string())?;
    // Finish every source admission and reduction before creating either file.
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, sample_count, options.block, options.full_scale_pa,
    )?;
    let parts: Vec<_> = sources.iter().zip(renderer.part_info()).map(|(source, aligned)| {
        let o = aligned.observation;
        let scope = if source.kind == Kind::Reed {
            ",\"observation_scope\":\"bore-pressure-plus-compact-jet-proxy; not exterior microphone\""
        } else if !source.exterior_json.is_empty() {
            source.exterior_json.as_str()
        } else if source.kind == Kind::Valve {
            ",\"observation_scope\":\"internal coupled tube pressure; not an exterior microphone\""
        } else { "" };
        format!("{{\"kind\":\"{}\",\"schema\":\"{}\",\"blake3\":\"{}\",\
            \"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"ratio\":{},\
            \"filter\":\"{}\",\"filter_delay_output_samples\":{},\"alignment_delay_samples\":{},\
            \"first_output_source_index\":{},\"source_full_scale_pa\":{:e},\"source_pcm_scale_applied\":false,\
            \"components\":{},\"modes\":{},\"compiled_controls\":{}{scope}}}",
            source.kind.label(), source.schema, source.hash.to_hex(), source.rate, source.samples,
            o.ratio, o.filter_profile, o.delay_output_samples, aligned.alignment_delay_samples,
            o.first_output_source_index, source.source_full_scale_pa, source.components, source.modes, source.controls)
    }).collect();
    let hash = rendered.hash.to_hex();
    let provenance = format!("{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"ensemble-input\",\
        \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{},\"full_scale_pa\":{:e},\
        \"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{hash}\",\
        \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
        \"ensemble\":{{\"sum_order\":\"input-order\",\"common_delay_output_samples\":{},\
        \"delay_compensated\":false,\"initial_history\":\"zero\",\"tail\":\"no-flush-declared-window\",\
        \"parts\":[{}],\"model_scope\":\"independent source pressures; common observer and time origin required; no new mechanical coupling or physical-validation claim\"}}}}",
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
