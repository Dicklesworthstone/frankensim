//! Explicit broadband modal file sub-lane using the common streaming PCM owner.
use super::super::{RATE, create_outputs, json_string, stream_output};
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::render::schedule::force::file::broadband::BroadbandModalPerformance;
use std::{io::Write, path::Path};

pub(super) fn run(
    bytes: &[u8],
    output: &Path,
    sidecar: &Path,
    block: usize,
    decimate: bool,
) -> Result<(), String> {
    let performance = BroadbandModalPerformance::from_bytes(bytes, block)?;
    let info = performance.info();
    if (!decimate && info.sample_rate_hz != RATE) || (decimate && info.sample_rate_hz <= RATE) {
        return Err("broadband mechanics require the same explicit --decimate clock opt-in as ordinary modal files".into());
    }
    let compiled_controls = performance.compiled_controls();
    let reports=performance.reports().iter().map(|r|format!(
        "{{\"producer_id\":{},\"authority\":\"EstimateOnly\",\"producer_relative_error_estimate\":{:e},\"fit_order\":{},\"fit_iterations\":{},\"delay_kernel_half_len\":{},\"delay_samples\":{:e},\"held_out_max_rel\":{:e},\"tolerance\":{:e},\"input_batches\":{},\"constant_narrow_band_compatibility\":{}}}",
        json_string(&r.source_id),r.producer_relative_error_estimate,r.controls.fit_order,
        r.controls.fit_iterations,r.controls.delay_kernel_half_len,r.delay_samples,
        r.held_out_max_rel_error,r.controls.max_held_out_rel_error,r.input_batches,r.constant_narrow_band_compatibility
    )).collect::<Vec<_>>().join(",");
    let mut renderer = DecimatedRenderer::new(performance, info.sample_rate_hz, RATE, block)
        .map_err(|e| e.to_string())?;
    let samples = usize::try_from(
        renderer
            .output_samples_for(info.samples)
            .map_err(|e| e.to_string())?,
    )
    .map_err(|_| "sample count exceeds host")?;
    let observation = renderer.info();
    // No artifact is created until complete parsing, fitting, held-out gating,
    // force compilation and output-clock admission have succeeded.
    let (mut audio, mut metadata) = create_outputs(output, sidecar)?;
    let rendered = stream_output::render_waveform(
        &mut renderer,
        &mut audio,
        samples,
        block,
        info.full_scale_pa,
    )?;
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"modal-broadband\",\"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\"full_scale_pa\":{:e},\"clipped_samples\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{}\",\"modal_input\":{{\"schema\":\"{}\",\"blake3\":\"{}\",\"voices\":{},\"modes\":{},\"force_events\":{},\"compiled_controls\":{compiled_controls},\"broadband\":[{reports}],\"integrator\":\"augmented continuous mode/filter exact ZOH\",\"delay_scope\":\"bandlimited interpolation of sampled pressure; zero acoustic history\",\"mechanics_sample_rate_hz\":{},\"observation_ratio\":{},\"observation_filter\":\"{}\",\"observation_delay_output_samples\":{},\"radiation\":\"producer transfer samples, gas-invariant\",\"model_scope\":\"authored modal model and declared computed receiver; EstimateOnly\"}}}}",
        info.full_scale_pa,
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        rendered.hash.to_hex(),
        info.schema,
        info.input_hash.to_hex(),
        info.voices,
        info.modes,
        info.force_events,
        info.sample_rate_hz,
        observation.ratio,
        observation.filter_profile,
        observation.delay_output_samples
    );
    writeln!(metadata, "{provenance}")
        .and_then(|_| metadata.flush())
        .map_err(|e| format!("broadband sidecar write failed: {e}"))?;
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"modal-broadband\",\"wav\":{},\"samples\":{samples},\"clipped\":{},\"peak_pa\":{:e},\"rms_pa\":{:e},\"wav_blake3\":\"{}\",\"input_blake3\":\"{}\",\"broadband\":[{reports}]}}",
        json_string(&output.display().to_string()),
        rendered.clipped,
        rendered.peak_pa,
        rendered.rms_pa,
        rendered.hash.to_hex(),
        info.input_hash.to_hex()
    );
    Ok(())
}
