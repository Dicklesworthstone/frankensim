//! Offline command output: one pressure block, one PCM block, then a bounded
//! readback of the finalized WAV for its existing replay identity. The reed's
//! baked characteristic-line storage is unchanged; only output staging is bounded.

use std::io::{BufWriter, Read, Seek, SeekFrom, Write};

use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::pcm_wav::stream::{MAX_PCM16_WAV_SAMPLES, Pcm16WavStream};
use fs_couple::render::schedule::ScheduledRenderer;

use super::{RATE, WAV_HASH_DOMAIN};

pub(super) struct RenderedWaveform {
    pub clipped: u64,
    pub peak_pa: f64,
    pub rms_pa: f64,
    pub hash: ContentHash,
}

pub(super) fn render_waveform<W: Read + Write + Seek>(
    renderer: &mut ScheduledRenderer,
    output: &mut W,
    samples: usize,
    block: usize,
    full_scale_pa: f64,
) -> Result<RenderedWaveform, String> {
    if samples == 0 || samples as u128 > u128::from(MAX_PCM16_WAV_SAMPLES)
        || block == 0 || block > renderer.context().max_block_len()
    {
        return Err("render output needs a positive RIFF-sized history and an admitted block".into());
    }
    renderer.validate_sample_rate(RATE).map_err(|e| e.to_string())?;
    if renderer.samples_rendered() != 0 {
        return Err("the command output starts at sample zero; use the stream API for continuation".into());
    }
    if output.stream_position().map_err(|e| e.to_string())? != 0 {
        return Err("the command WAV must start at byte zero".into());
    }
    // Buffer tiny --block writes without changing the physics partition. The
    // stream owns the canonical PCM scaling and the single finalized header.
    let buffered = BufWriter::new(&mut *output);
    let mut stream = Pcm16WavStream::new(buffered, RATE, full_scale_pa, block)
        .map_err(|e| e.to_string())?;
    let mut pressure = vec![0.0_f64; block.min(samples)];
    let mut peak = 0.0_f64;
    let mut squares = 0.0_f64;
    let mut completed = 0;
    while completed < samples {
        let count = pressure.len().min(samples - completed);
        let chunk = &mut pressure[..count];
        renderer.block(chunk).map_err(|e| format!("render refused: {e}"))?;
        // Retain sample order, not per-block partial sums: callback boundaries
        // must not change the statistics used in provenance.
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
    let (buffered, summary) = stream.finish().map_err(|e| e.to_string())?;
    // finish flushed all payload/header bytes. Release the exclusive file
    // borrow before reading the SAME handle; no path-based reopen can substitute
    // another file between production and hashing.
    drop(buffered);
    output.seek(SeekFrom::Start(0)).map_err(|e| format!("WAV hash seek failed: {e}"))?;
    let mut remaining = 44 + 2 * summary.samples;
    let mut bytes = [0_u8; 8192];
    let mut hasher = DomainHasher::new(WAV_HASH_DOMAIN);
    while remaining > 0 {
        let count = remaining.min(bytes.len() as u64) as usize;
        output.read_exact(&mut bytes[..count]).map_err(|e| format!("WAV hash read failed: {e}"))?;
        hasher.update(&bytes[..count]);
        remaining -= count as u64;
    }
    Ok(RenderedWaveform {
        clipped: summary.clipped_samples,
        peak_pa: peak,
        rms_pa: (squares / samples as f64).sqrt(),
        hash: hasher.finalize(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
    use fs_couple::pcm_wav::encode_pcm16_wav;
    use fs_couple::render::{ControlDelta, ModalStringVoice, RenderContext, RenderVoice};
    use fs_couple::render::schedule::ScheduledControl;

    fn image(rate: u32) -> ModalAcousticTimeModel {
        ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
            angular_frequency_rad_s: core::f64::consts::TAU * 440.0,
            damping_ratio: 0.02,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        }], ModalAcousticTimeBudget::audible_reference()).unwrap()
    }

    fn renderer(rate: u32) -> ScheduledRenderer {
        let voice = ModalStringVoice::new(image(rate), vec![1.0]).unwrap();
        ScheduledRenderer::new(RenderContext::new(vec![RenderVoice::ModalString(voice)], 512),
            vec![ScheduledControl { sample: 37, delta: ControlDelta::SetModalForce {
                voice: 0, mode: 0, force_n_per_sqrt_kg: 0.0,
            }}], 1).unwrap()
    }

    #[test]
    fn command_output_matches_owned_encoding_and_identity_with_partition_independent_statistics() {
        let mut direct = image(RATE);
        let pressure: Vec<_> = (0..257).map(|i| direct.step(&[if i < 37 { 1.0 } else { 0.0 }])
            .unwrap().observer_pressure_pa).collect();
        let (expected, clips) = encode_pcm16_wav(&pressure, RATE, 0.0001).unwrap();
        let expected_hash = fs_blake3::hash_domain(WAV_HASH_DOMAIN, &expected);
        let expected_peak = pressure.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let expected_rms = (pressure.iter().map(|v| v*v).sum::<f64>() / 257.0).sqrt();
        let mut previous = None;
        for block in [1, 7, 37, 64, 512] {
            let mut output = Cursor::new(Vec::new());
            let report = render_waveform(&mut renderer(RATE), &mut output, 257, block, 0.0001).unwrap();
            assert_eq!(output.into_inner(), expected);
            assert_eq!(report.hash, expected_hash);
            assert_eq!(report.clipped, clips as u64);
            assert_eq!(report.peak_pa, expected_peak);
            assert!((report.rms_pa - expected_rms).abs() <= 1e-14 * expected_rms);
            let bits = (report.peak_pa.to_bits(), report.rms_pa.to_bits());
            if let Some(previous) = previous { assert_eq!(bits, previous); }
            previous = Some(bits);
        }
    }

    #[test]
    fn command_admission_refuses_without_creating_payload_or_advancing_models() {
        for (samples, block, rate) in [(0, 32, RATE), (10, 0, RATE), (10, 513, RATE), (10, 32, 44_100)] {
            let mut renderer = renderer(rate);
            let mut output = Cursor::new(Vec::new());
            assert!(render_waveform(&mut renderer, &mut output, samples, block, 1.0).is_err());
            assert_eq!(renderer.samples_rendered(), 0);
            assert!(output.into_inner().is_empty());
        }
    }

    struct Unreadable(Cursor<Vec<u8>>);
    impl Read for Unreadable {
        fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("injected readback failure"))
        }
    }
    impl Write for Unreadable {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> { self.0.write(bytes) }
        fn flush(&mut self) -> std::io::Result<()> { self.0.flush() }
    }
    impl Seek for Unreadable {
        fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> { self.0.seek(position) }
    }

    #[test]
    fn unreadable_output_cannot_publish_a_replay_hash() {
        let mut output = Unreadable(Cursor::new(Vec::new()));
        let result = render_waveform(&mut renderer(RATE), &mut output, 257, 37, 1.0);
        assert!(result.err().unwrap().contains("WAV hash read failed"));
    }
}
