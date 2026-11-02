//! Continuous physical tone forces on one original v2 coupled assembly.
use super::super::{RATE, create_outputs, json_string, stream_output};
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::render::schedule::force::file::coupled_tones::CoupledTonePerformance;
use std::{io::Write, path::Path};

pub(super) fn run(
    bytes: &[u8],
    output: &Path,
    sidecar: &Path,
    block: usize,
    decimate: bool,
) -> Result<(), String> {
    let performance = CoupledTonePerformance::from_bytes(bytes, block)?;
    let info = performance.info();
    if (!decimate && info.sample_rate_hz != RATE) || (decimate && info.sample_rate_hz <= RATE) {
        return Err("coupled tone mechanics require the same explicit --decimate clock opt-in as ordinary modal files".into());
    }
    let controls = performance.compiled_controls();
    let tones = performance.continuous_force_tones();
    let mut renderer = DecimatedRenderer::new(performance, info.sample_rate_hz, RATE, block)
        .map_err(|e| e.to_string())?;
    let samples = usize::try_from(
        renderer
            .output_samples_for(info.samples)
            .map_err(|e| e.to_string())?,
    )
    .map_err(|_| "sample count exceeds host")?;
    let observation = renderer.info();
    // Parse, project, preload, prepare all continuous work/heat and admit the
    // observation clock before the shared PCM owner creates any output file.
    let (mut audio, mut metadata) = create_outputs(output, sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer,
        &mut audio,
        samples,
        block,
        info.full_scale_pa,
    )?;
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"modal-coupled-tones\",\"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{}\",\"modal_input\":{{\"schema\":\"{}\",\"blake3\":\"{}\",\"voices\":{},\"modes\":{},\"connections\":{},\"force_events\":{},\"compiled_controls\":{controls},\"continuous_force_tones\":{tones},\"force_phase_origin_s\":0,\"integrator\":\"complete coupled mechanics and continuous forcing exact exponential\",\"heat\":\"independent original loss channels including cross-tone work\",\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"observation_ratio\":{},\"observation_filter\":\"{}\",\"observation_delay_output_samples\":{},\"radiation\":\"authored real instantaneous pressure per modal velocity; gas-invariant\",\"model_scope\":\"authored linear reduced assembly; no physical-validation claim\"}}}}",
        info.full_scale_pa,
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        rendered.hash.to_hex(),
        info.schema,
        info.input_hash.to_hex(),
        info.voices,
        info.modes,
        info.connections,
        info.force_events,
        info.sample_rate_hz,
        info.samples,
        observation.ratio,
        observation.filter_profile,
        observation.delay_output_samples,
    );
    writeln!(metadata, "{provenance}")
        .and_then(|_| metadata.flush())
        .map_err(|e| format!("coupled tone sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"modal-coupled-tones\",\"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{}\",\"input_blake3\":\"{}\"}}",
        json_string(&output.display().to_string()),
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        rendered.hash.to_hex(),
        info.input_hash.to_hex(),
    );
    Ok(())
}
