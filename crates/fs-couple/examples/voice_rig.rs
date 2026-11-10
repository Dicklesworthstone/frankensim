//! Research render of one voice performance card at the clock the card declares, for flow-rig comparisons whose
//! printed boundary fixes the clock: `termination unflanged_thin_wall` is refused by fs-duct above the lip section's
//! first cut-on (`ka` 3.83), so a rig with a wide lip runs below the site's 48 kHz. `music_render voice` keeps its
//! 48 kHz site clock and is not changed; nothing written here is a site clip. The card is parsed and built by the same
//! [`VoiceDescription`] owner, and the samples stream through the same [`Pcm16WavStream`] in the same 512-sample blocks,
//! so a 48 kHz card gives the bytes `music_render voice` writes for it.
//!
//! usage: `cargo run --release -p fs-couple --example voice_rig -- INPUT.performance OUT.wav`
//! The output must not exist. One JSON line goes to stdout: the clock, the sample count, clipping, the pressure peak and
//! RMS, the glottal flow extremes and, for a card with `glottal_orifice`, the mean glottal flow over the last 2 s.
//! A card with an airway below the glottis logs that line's length, rung count and passivity gain to stderr first.

use std::fs::File;
use std::io::{BufWriter, Read};

use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::stream::{MAX_PCM16_WAV_SAMPLES, Pcm16WavStream};
use fs_couple::voice_performance::{MAX_VOICE_PERFORMANCE_BYTES, VoiceDescription};

const BLOCK: usize = 512;

fn main() {
    if let Err(why) = run() {
        eprintln!("{{\"suite\":\"voice-rig\",\"verdict\":\"refused\",\"what\":{why:?}}}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input, output] = args.as_slice() else {
        return Err("usage: voice_rig INPUT.performance OUT.wav".into());
    };
    let mut bytes = Vec::new();
    File::open(input)
        .map_err(|e| format!("voice input open failed: {e}"))?
        .take((MAX_VOICE_PERFORMANCE_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("voice input read failed: {e}"))?;
    // No gas registry: a card that needs one is refused here, as it is by the bin without `--gas-component`.
    let description = VoiceDescription::from_bytes_with_gas_sources(&bytes, None)?;
    let rate = description.output_rate_hz;
    let samples = usize::try_from(description.samples).map_err(|_| "sample count exceeds this host".to_string())?;
    if samples == 0 || samples as u128 > u128::from(MAX_PCM16_WAV_SAMPLES) {
        return Err("the render needs a positive RIFF-sized history".into());
    }
    let mut performance = description.build()?;
    if let Some(airway) = &performance.info().subglottal_airway {
        let line = &airway.reflectance;
        let gain = line.passivity_gain.map_or_else(
            || "null".to_string(),
            |g| format!("{{\"passes\":{},\"top_band_db\":{:e}}}", g.passes, g.top_band_db),
        );
        eprintln!(
            "{{\"suite\":\"voice-rig\",\"airway_line_n_fft\":{},\"airway_line_rungs\":{},\"airway_line_passivity_gain\":{gain}}}",
            line.n_fft,
            line.rungs_tried.len()
        );
    }
    performance.validate_sample_rate(rate).map_err(|e| e.to_string())?;
    performance.validate_sample_count(samples as u64).map_err(|e| e.to_string())?;
    if BLOCK > performance.max_block_len() {
        return Err("the voice performance does not admit a 512-sample block".into());
    }
    let file = File::create_new(output).map_err(|e| format!("output create failed (it must not exist): {e}"))?;
    let mut stream = Pcm16WavStream::new(BufWriter::new(file), rate, description.full_scale_pa, BLOCK)
        .map_err(|e| e.to_string())?;
    let mut pressure = vec![0.0_f64; BLOCK.min(samples)];
    let (mut peak, mut squares, mut completed) = (0.0_f64, 0.0_f64, 0);
    while completed < samples {
        let count = pressure.len().min(samples - completed);
        let chunk = &mut pressure[..count];
        performance.block(chunk).map_err(|e| format!("render refused: {e}"))?;
        for &value in chunk.iter() {
            peak = peak.max(value.abs());
            squares += value * value;
        }
        if !peak.is_finite() || !squares.is_finite() {
            return Err("pressure statistics exceeded the finite range; output is incomplete".into());
        }
        stream.write_block(chunk).map_err(|e| e.to_string())?;
        completed += count;
    }
    let (_, summary) = stream.finish().map_err(|e| e.to_string())?;
    let stats = performance.stats();
    let rig = stats.rig_mean_glottal_flow_m3_s.map_or_else(|| "null".to_string(), |u| format!("{u:e}"));
    println!(
        "{{\"suite\":\"voice-rig\",\"output_rate_hz\":{rate},\"site_clip\":false,\"samples\":{samples},\
         \"clipped\":{},\"peak_pa\":{peak:e},\"rms_pa\":{:e},\"peak_glottal_flow_m3_s\":{:e},\
         \"min_glottal_flow_m3_s\":{:e},\"mean_glottal_flow_last_2_s_m3_s\":{rig}}}",
        summary.clipped_samples,
        (squares / samples as f64).sqrt(),
        stats.peak_glottal_flow_m3_s,
        stats.min_glottal_flow_m3_s,
    );
    Ok(())
}
