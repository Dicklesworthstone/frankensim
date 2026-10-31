//! Glottal-source plus tract-area-function command; reuses the streaming output owner.

use super::{RATE, create_outputs, json_string, stream_output};
use fs_blake3::DomainHasher;
use fs_couple::voice_performance::{
    MAX_VOICE_PERFORMANCE_BYTES, VOICE_PERFORMANCE_SCHEMA, VoiceDescription,
};
use fs_material::gas_species::input::ThermochemicalSources;
use std::io::{Read, Write};
use std::path::Path;

const INPUT_HASH_DOMAIN: &str = "org.frankensim.fs-couple.voice-performance.v1";

fn options(args: &[String]) -> Result<(&str, &str, usize), String> {
    let mut paths = Vec::new();
    let mut block = None;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--block" => {
                if block.is_some() { return Err("--block may only be supplied once".into()); }
                block = Some(iter.next().and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| (1..=65_536).contains(v))
                    .ok_or_else(|| "--block needs an integer in 1..=65536".to_string())?);
            }
            value if value.starts_with('-') => return Err(format!(
                "unsupported voice option {value:?}; the fold, tract, gas, duration, pressure and microphone belong to the file"
            )),
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render voice INPUT.performance OUT.wav [--block N]".into());
    };
    Ok((*input, *output, block.unwrap_or(512)))
}

pub(super) fn run_with_gas_sources(
    args: &[String],
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(), String> {
    let (input, output, block) = options(args)?;
    let output = Path::new(output);
    let sidecar = output.with_extension("provenance.json");
    if output == sidecar.as_path() || output.exists() || sidecar.exists() {
        return Err("output/sidecar paths must be distinct and must not already exist".into());
    }
    let file = std::fs::File::open(input).map_err(|e| format!("voice input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_VOICE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("voice input read failed: {e}"))?;
    let description = VoiceDescription::from_bytes_with_gas_sources(&bytes, sources)?;
    if description.output_rate_hz != RATE {
        return Err("music_render requires an explicitly declared 48000 Hz output clock; no implicit resampling".into());
    }
    let samples = usize::try_from(description.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let mut performance = description.build()?;
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut performance, &mut audio, samples, block, description.full_scale_pa,
    )?;
    let mut hasher = DomainHasher::new(INPUT_HASH_DOMAIN);
    hasher.update(&bytes);
    let input_hash = hasher.finalize().to_hex();
    let wav_hash = rendered.hash.to_hex();
    let info = performance.info();
    let stats = performance.stats();
    let mic = description.microphone.position_m;
    let gas_sources =
        super::thermochemical::provenance(core::slice::from_ref(description.gas_record()));
    // Rigid-end limit of this tract, mouth load and wall; see VoiceInfo::passive_formants_hz.
    let passive = match &info.passive_formants_hz {
        Ok(formants) => format!(
            "[{}]",
            formants.iter().map(|f| format!("{f:.1}")).collect::<Vec<_>>().join(",")
        ),
        Err(why) => format!("{{\"refused\":{}}}", json_string(why)),
    };
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"voice-input\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
         \"voice_input\":{{\"schema\":\"{VOICE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\
         \"mechanics_rate_hz\":{},\"two_mass\":{},\"sections\":{},\"tract_length_m\":{:e},\
         \"glottal_radius_m\":{:e},\"mouth_radius_m\":{:e},\"sound_speed_m_s\":{:e},\
         \"yielding_wall\":{},\"tract_source\":{},\"tract_licence\":{},\
         \"microphone_m\":[{:e},{:e},{:e}],\"mouth_flow_taps\":{},\"mouth_flow_rungs\":{},\
         \"passive_tract_formants_hz\":{passive},\
         \"peak_glottal_flow_m3_s\":{:e},\"closed_fraction\":{:e},\
         \"max_force_passes\":{},\"mean_force_passes\":{:e},\
         \"model_scope\":\"lumped fold pHS on a plane-wave series tract with TMM viscothermal loss and a declared radiation load; mouth flow radiated through an infinite-baffle circular opening, no head diffraction, nasal branch or side cavity; the supraglottal pressure is never the output\"}}{gas_sources}}}",
        description.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        info.mechanics_rate_hz, info.two_mass, info.sections, info.tract_length_m,
        info.glottal_radius_m, info.mouth_radius_m, info.sound_speed_m_s,
        info.has_wall, json_string(&description.tract_source), json_string(&description.tract_licence),
        mic[0], mic[1], mic[2], info.mouth_flow_taps, info.mouth_flow_rungs,
        stats.peak_glottal_flow_m3_s, stats.closed_fraction,
        stats.max_force_passes, stats.mean_force_passes,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"voice-input\",\
         \"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"peak_glottal_flow_m3_s\":{:e},\"closed_fraction\":{:e},\
         \"max_force_passes\":{},\"mean_force_passes\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, rendered.peak_pa,
        rendered.rms_pa, stats.peak_glottal_flow_m3_s, stats.closed_fraction,
        stats.max_force_passes, stats.mean_force_passes,
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_string()).collect() }
    #[test]
    fn voice_file_owns_physics_and_clock_while_block_size_remains_host_owned() {
        for flag in ["--seconds", "--schedule", "--full-scale-pa", "--frequency", "--vowel"] {
            assert!(options(&args(&["input", "output", flag, "1"])).is_err());
        }
        for value in ["0", "65537", "-1", "NaN"] {
            assert!(options(&args(&["input", "output", "--block", value])).is_err());
        }
        assert!(options(&args(&["input", "output", "--block", "37", "--block", "64"])).is_err());
        assert!(options(&args(&["input"])).is_err());
        let a = args(&["input", "output", "--block", "37"]);
        assert_eq!(options(&a).unwrap(), ("input", "output", 37));
    }
}
