//! Authored reed/duct pressure to WAV; all mechanics and PCM stay in their owners.
use std::io::{Read, Write};
use std::path::Path;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::render::schedule::reed::{ReedPerformance, MAX_REED_PERFORMANCE_BYTES, REED_PERFORMANCE_SCHEMA};
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

pub(super) fn run(args: &[String]) -> Result<(), String> {
    let (input, output, block, decimate) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("wind input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_REED_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("wind input read failed: {e}"))?;
    let performance = ReedPerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
    let info = performance.info();
    if !decimate && info.sample_rate_hz != RATE {
        return Err("non-48000-Hz mechanics require explicit --decimate; no implicit resampling".into());
    }
    if decimate && info.sample_rate_hz <= RATE {
        return Err("--decimate requires mechanics above 48000 Hz; ordinary 48000-Hz input needs no conversion".into());
    }
    // Keep the finite wrapper inside the observer; extracting a raw scheduler
    // would discard its end-of-performance admission. The ratio-one path is
    // arithmetic identity, while higher rates filter ONLY the observed pressure.
    let mut renderer = DecimatedRenderer::new(performance, info.sample_rate_hz, RATE, block)
        .map_err(|e| e.to_string())?;
    let samples = renderer.output_samples_for(info.samples).map_err(|e| e.to_string())?;
    if samples == 0 { return Err("wind source needs at least one complete output interval".into()); }
    renderer.validate_sample_count(samples).map_err(|e| e.to_string())?;
    let samples = usize::try_from(samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let observation = renderer.info();
    // Complete source/clock/window admission before creating either artifact.
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, samples, block, info.full_scale_pa,
    )?;
    let input_hash = info.input_hash.to_hex();
    let wav_hash = rendered.hash.to_hex();
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"reed-input\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
         \"reed_input\":{{\"schema\":\"{REED_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\
         \"segments\":{},\"tone_holes\":{},\"massive_reed\":{},\"gesture_events\":{},\"compiled_controls\":{},\
         \"model_scope\":\"authored reed primitives and static duct; bore-pressure plus compact-jet proxy; not a calibrated exterior microphone or material-identification claim\"}},\
         \"observation\":{{\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\
         \"output_sample_rate_hz\":{RATE},\"ratio\":{},\"filter\":\"{}\",\
         \"delay_output_samples\":{},\"first_output_source_index\":{},\
         \"initial_history\":\"zero\",\"delay_compensated\":false,\"tail\":\"no-flush-declared-window\"}}}}",
        info.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        info.segments, info.tone_holes, info.massive_reed, info.gesture_events, info.compiled_controls,
        info.sample_rate_hz, info.samples, observation.ratio, observation.filter_profile,
        observation.delay_output_samples, observation.first_output_source_index,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!("{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"reed-input\",\
        \"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
        \"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa, rendered.rms_pa);
    Ok(())
}
