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
    let envelopes = performance.envelope_count();
    let sweeps = performance.frequency_sweep_count();
    let force_error = performance.maximum_sweep_force_error_bound_n();
    let reconstruction = json_string(if sweeps == 0 {
        "exact prescribed constant or affine carrier coefficients"
    } else {
        "affine carrier reconstruction of integrated source-frequency curves; force bound per original tone"
    });
    let (broadband, radiation) = match performance.broadband_report() {
        None => (
            "null".to_owned(),
            "authored real instantaneous pressure per modal velocity; gas-invariant",
        ),
        Some(r) => (
            format!(
                "{{\"producer_id\":{},\"authority\":\"EstimateOnly\",\"producer_relative_error_estimate\":{:e},\"fit_order\":{},\"fit_iterations\":{},\"delay_kernel_half_len\":{},\"delay_samples\":{:e},\"held_out_max_rel\":{:e},\"tolerance\":{:e},\"input_batches\":{},\"constant_narrow_band_compatibility\":{},\"delay_scope\":\"bandlimited interpolation of sampled pressure; zero receiver history\"}}",
                json_string(&r.source_id),
                r.producer_relative_error_estimate,
                r.controls.fit_order,
                r.controls.fit_iterations,
                r.controls.delay_kernel_half_len,
                r.delay_samples,
                r.held_out_max_rel_error,
                r.controls.max_held_out_rel_error,
                r.input_batches,
                r.constant_narrow_band_compatibility,
            ),
            "complete coupled trajectory through sampled causal pressure transfer; gas-invariant",
        ),
    };
    let radiation = json_string(radiation);
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
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"modal-coupled-tones\",\"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{}\",\"modal_input\":{{\"schema\":\"{}\",\"blake3\":\"{}\",\"voices\":{},\"modes\":{},\"connections\":{},\"force_events\":{},\"compiled_controls\":{controls},\"continuous_force_tones\":{tones},\"physical_force_envelopes\":{envelopes},\"physical_force_frequency_sweeps\":{sweeps},\"maximum_tone_force_reconstruction_error_n\":{force_error:e},\"source_reconstruction\":{reconstruction},\"force_phase_origin_s\":0,\"integrator\":\"complete coupled mechanics and reconstructed continuous forcing exact exponential\",\"heat\":\"independent original loss channels including cross-tone work\",\"mechanics_sample_rate_hz\":{},\"mechanics_samples\":{},\"observation_ratio\":{},\"observation_filter\":\"{}\",\"observation_delay_output_samples\":{},\"broadband\":{broadband},\"radiation\":{radiation},\"model_scope\":\"authored linear reduced assembly; no physical-validation claim\"}}}}",
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
