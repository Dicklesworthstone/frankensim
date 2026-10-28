//! Geometry-driven plate command; reuse the existing streaming output owner.

use std::io::{Read, Write};
use std::path::Path;
use fs_material::gas_species::input::ThermochemicalSources;
use fs_matdb::NormalizedMaterialCardPack;
use fs_couple::render::plate::file::{
    MAX_PLATE_PERFORMANCE_BYTES, PLATE_PERFORMANCE_SCHEMA, PlatePerformance,
};
use super::{RATE, create_outputs, json_string, stream_output};

#[path = "material_packs.rs"]
mod material_packs;
pub(super) fn section_pack_loader(parent: &Path)
    -> impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String> + '_
{
    material_packs::section_pack_loader(parent)
}

/// `,"section_sources":[...]`: each section's pack, hash, state and point, or its fixture label.
pub(super) fn section_sources_json(sources: &[String]) -> String {
    format!(",\"section_sources\":[{}]", sources.iter().map(|s| format!("\"{s}\"")).collect::<Vec<_>>().join(","))
}

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
                "unsupported plate option {value:?}; geometry, duration, force and full-scale belong to the file"
            )),
            value => paths.push(value),
        }
    }
    let [input, output] = paths.as_slice() else {
        return Err("usage: music_render plate INPUT.performance OUT.wav [--block N]".into());
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
    let file = std::fs::File::open(input).map_err(|e| format!("plate input open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_PLATE_PERFORMANCE_BYTES + 1) as u64).read_to_end(&mut bytes)
        .map_err(|e| format!("plate input read failed: {e}"))?;
    let parent = Path::new(input).parent().unwrap_or(Path::new("."));
    let performance = PlatePerformance::from_bytes_with_all_sources(&bytes, block, sources, section_pack_loader(parent))
        .map_err(|e| e.to_string())?;
    let section_sources = section_sources_json(performance.section_provenance());
    let gas_sources = super::thermochemical::provenance(core::slice::from_ref(performance.gas_record()));
    let info = performance.info();
    if info.sample_rate_hz != RATE {
        return Err("music_render requires an explicitly declared 48000 Hz input; no implicit resampling".into());
    }
    let samples = usize::try_from(info.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    let mut renderer = performance.into_renderer();
    let (mut audio, mut metadata) = create_outputs(output, &sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer, &mut audio, samples, block, info.full_scale_pa,
    )?;
    let input_hash = info.input_hash.to_hex();
    let wav_hash = rendered.hash.to_hex();
    let radiation_load = if info.baffled_piston_loading { "baffled-piston" } else { "none" };
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"plate-input\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\
         \"wav_blake3\":\"{wav_hash}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\",\
         \"plate_input\":{{\"schema\":\"{PLATE_PERFORMANCE_SCHEMA}\",\"blake3\":\"{input_hash}\",\
         \"nodes\":{},\"triangles\":{},\"sections\":{},\"requested_modes\":{},\
         \"retained_modes\":{},\"force_events\":{},\
         \"ambient_gas\":{{\"temperature_k\":{:e},\"pressure_pa\":{:e},\"density_kg_m3\":{:e},\
         \"sound_speed_m_s\":{:e},\"water_mole_fraction\":{:e}}},\
         \"listener_m\":{:e},\"propagation_delay_s\":{:e},\"radiation_load\":\"{radiation_load}\"{section_sources},\
         \"model_scope\":\"linear flat DKT plate; sections from their named packs or fixture labels; compact baffled observer retarded by r/c; optional per-mode piston reaction uses the same gas and clock; no mutual-radiation, continuum-convergence or physical-validation claim\"}}{gas_sources}}}",
        info.full_scale_pa, rendered.clipped, rendered.peak_pa, rendered.rms_pa,
        info.nodes, info.triangles, info.sections, info.requested_modes, info.retained_modes,
        info.force_events, info.gas.temperature, info.gas.pressure, info.gas.density,
        info.gas.sound_speed, info.gas.water_mole_fraction, info.listener_m, info.propagation_delay_s,
    );
    writeln!(metadata, "{provenance}").and_then(|_| metadata.flush())
        .map_err(|e| format!("sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"plate-input\",\
         \"wav\":{},\"samples\":{samples},\"clipped\":{},\"retained_modes\":{},\
         \"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{wav_hash}\",\"input_blake3\":\"{input_hash}\"}}",
        json_string(&output.display().to_string()), rendered.clipped, info.retained_modes,
        rendered.peak_pa, rendered.rms_pa,
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|s| (*s).to_string()).collect() }
    #[test]
    fn plate_file_owns_physics_and_clock_while_block_size_remains_host_owned() {
        for flag in ["--seconds", "--schedule", "--full-scale-pa", "--frequency"] {
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
