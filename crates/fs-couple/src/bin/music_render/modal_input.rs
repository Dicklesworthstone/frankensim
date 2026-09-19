//! File-driven reduced models, sharing the existing command's streaming output.

use std::io::{Read, Write};
use std::path::Path;
use fs_couple::render::schedule::force::file::{
    MAX_MODAL_PERFORMANCE_BYTES, MODAL_PERFORMANCE_SCHEMA, MODAL_CONTACT_PERFORMANCE_SCHEMA,
    MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA, ModalPerformance,
};
use super::{RATE, create_outputs, json_string, stream_output};

fn options(args: &[String]) -> Result<(&str, &str, usize), String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                let value = iter.next().and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?;
                block = Some(value);
            }
            value if value.starts_with('-') => return Err(format!(
                "unsupported modal option {value:?}; duration, full-scale and forces belong to the input file"
            )),
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render modal INPUT.performance OUT.wav [--block N]".into());
    };
    Ok((*input, *output, block.unwrap_or(512)))
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    let (input, output, block) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let input = std::fs::File::open(input).map_err(|e| format!("modal input open failed: {e}"))?;
    let mut bytes = Vec::new();
    input.take((MAX_MODAL_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("modal input read failed: {e}"))?;
    let performance = ModalPerformance::from_bytes(&bytes, block).map_err(|e| e.to_string())?;
    let info = performance.info();
    if info.sample_rate_hz != RATE {
        return Err("music_render requires an explicitly declared 48000 Hz input; no implicit resampling".into());
    }
    let samples = usize::try_from(info.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let mut renderer = performance.into_renderer();
    let compiled_controls = renderer.pending_controls().len();
    // Finish ALL input/model/clock admission before creating either artifact.
    // The existing output owner preserves stream scaling, short tails and hashes.
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, samples, block, info.full_scale_pa,
    )?;
    let input_hash = info.input_hash.to_hex();
    let wav_hash = rendered.hash.to_hex();
    let input_schema = info.schema;
    let mut coupling_provenance = if input_schema == MODAL_PERFORMANCE_SCHEMA { String::new() } else {
        format!(",\"connections\":{},\"coupling\":\"implicit-bilateral-spring-damper\"", info.connections)
    };
    if input_schema == MODAL_CONTACT_PERFORMANCE_SCHEMA {
        coupling_provenance.push_str(",\"contacts\":1,\"contact\":\"implicit-nonadhesive-power-law\"");
    } else if input_schema == MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA {
        coupling_provenance.push_str(&format!(
            ",\"contacts\":{},\"contact\":\"simultaneous-implicit-nonadhesive-power-law\"", info.contacts
        ));
    }
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"modal-input\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
         \"modal_input\":{{\"schema\":\"{input_schema}\",\"blake3\":\"{input_hash}\",\
         \"voices\":{},\"modes\":{},\"force_events\":{},\"compiled_controls\":{compiled_controls}{coupling_provenance},\
         \"model_scope\":\"authored reduced model; no physical-validation claim\"}}}}",
        info.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        info.voices, info.modes, info.force_events,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"modal-input\",\
         \"wav\":{},\"samples\":{samples},\"clipped\":{},\
         \"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa, rendered.rms_pa,
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_string()).collect() }

    #[test]
    fn model_file_owns_duration_scale_and_schedule_without_ignored_overrides() {
        for flag in ["--seconds", "--full-scale-pa", "--schedule", "--unknown"] {
            assert!(options(&args(&["input", "output", flag, "1"])).is_err());
        }
        for value in ["0", "65537", "-1", "NaN", "18446744073709551616"] {
            assert!(options(&args(&["input", "output", "--block", value])).is_err());
        }
        assert!(options(&args(&["input", "output", "--block"])).is_err());
        assert!(options(&args(&["input", "output", "--block", "37", "--block", "64"])).is_err());
        assert!(options(&args(&["input"])).is_err());
        assert!(options(&args(&["input", "output", "extra"])).is_err());
    }

    #[test]
    fn callback_size_is_the_only_optional_override() {
        let input = args(&["in.performance", "out.wav"]);
        assert_eq!(options(&input).unwrap(), ("in.performance", "out.wav", 512));
        let input = args(&["--block", "37", "in.performance", "out.wav"]);
        assert_eq!(options(&input).unwrap(), ("in.performance", "out.wav", 37));
    }
}
