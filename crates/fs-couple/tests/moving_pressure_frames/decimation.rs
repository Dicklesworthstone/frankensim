//! Real moving receivers retain native filter phase, channels and PCM bytes.

use super::*;
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::observation::{DecimatedFrameRenderer, DecimatedRenderer};
use std::io::{Seek, SeekFrom, Write};

fn array<'a>(
    left: &'a PreparedMovingReceiverSource<'a>,
    right: &'a PreparedMovingReceiverSource<'a>,
    samples: &'a [TimedCompactSource],
    rate: u32,
    source_frames: u64,
    capacity: usize,
) -> Frames<'a> {
    let cfg = MovingPressureConfig { sample_rate_hz: rate, samples: source_frames,
        max_block: capacity, ..config() };
    MovingPressureFrameRenderer::new(["approaching", "receding"],
        [renderer(left, samples, cfg).unwrap(), renderer(right, samples, cfg).unwrap()],
        budget(), &mut || false).unwrap()
}

#[test]
fn every_integer_ratio_matches_the_original_mono_decimators_bit_for_bit() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    const FRAMES: usize = 191;
    // Past the original filter startup support at every ratio. Each underlying
    // moving microphone has a distinct retarded source clock and pressure.
    for ratio in 1..=16 {
        let rate = 2000*ratio;
        let source_frames = FRAMES as u64 * u64::from(ratio);
        let cfg = MovingPressureConfig { sample_rate_hz: rate, samples: source_frames,
            max_block: 32, ..config() };
        let mut mono = [
            DecimatedRenderer::new(renderer(&left, &samples, cfg).unwrap(), rate, 2000, 23).unwrap(),
            DecimatedRenderer::new(renderer(&right, &samples, cfg).unwrap(), rate, 2000, 23).unwrap(),
        ];
        let mut expected = vec![0.0; FRAMES*2];
        for (channel, stream) in mono.iter_mut().enumerate() {
            let mut values = vec![0.0; FRAMES];
            for part in values.chunks_mut(23) { stream.block(part).unwrap(); }
            for (i, value) in values.into_iter().enumerate() { expected[2*i+channel] = value; }
        }
        for capacity in [1, 3, 32] {
            let block = if ratio == 1 { capacity.min(23) } else { 23 };
            let mut stream = DecimatedFrameRenderer::new(
                array(&left, &right, &samples, rate, source_frames, capacity), rate, 2000, block,
            ).unwrap();
            assert_eq!(stream.output_frames_for(source_frames).unwrap(), FRAMES as u64);
            assert_eq!(stream.info(), mono[0].info());
            assert_eq!(stream.info().first_output_source_index, ratio as usize-1);
            let mut actual = vec![0.0; FRAMES*2];
            for part in actual.chunks_mut(2*block) { stream.frame_block(part).unwrap(); }
            bits(&actual, &expected);
            assert_eq!(stream.frames_rendered(), FRAMES as u64);
            assert_eq!(stream.source().frames_rendered(), source_frames);
            assert!(stream.source().channels().iter().all(|r| r.samples_rendered()==source_frames));
        }
    }
}

#[test]
fn native_filtered_array_pcm_is_identical_across_callback_sizes() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    const FRAMES: usize = 181;
    let make = || DecimatedFrameRenderer::new(
        array(&left, &right, &samples, 12_000, 3*FRAMES as u64, 2), 12_000, 4000, 17,
    ).unwrap();
    let mut reference = make();
    let mut pressure = [0.0; 2*FRAMES];
    for part in pressure.chunks_mut(34) { reference.frame_block(part).unwrap(); }
    let (expected, clips) = encode_pcm16_wav_interleaved(&pressure, 4000, 2, 0.02).unwrap();
    for block in [1, 7, 17] {
        let mut renderer = make();
        let mut stream = Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 4000, 2, 0.02, block).unwrap();
        assert_eq!(render_pressure_frames_pcm16(&mut renderer, &mut stream,
            &CancelGate::new_clock_free(), &mut vec![0.0; 2*block], FRAMES as u64).unwrap(),
            WavFrameProgress::Completed { frames: FRAMES as u64 });
        let (sink, summary) = stream.finish().unwrap();
        assert_eq!(sink.into_inner(), expected);
        assert_eq!((summary.channels, summary.frames, summary.samples, summary.clipped_samples),
            (2, FRAMES as u64, 2*FRAMES as u64, clips as u64));
    }
}

#[test]
fn bad_groups_rates_and_missing_prehistory_refuse_before_advancing_motion() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let fresh = || array(&left, &right, &samples, 6000, 101, 3);
    for (input, output, block) in [(6000, 0, 5), (6000, 4000, 5), (6000, 12_000, 5),
        (6000, 300, 5), (6001, 2000, 5), (6000, 2000, 0), (6000, 6000, 4)] {
        assert!(DecimatedFrameRenderer::new(fresh(), input, output, block).is_err());
    }
    let mut started = fresh();
    started.frame_block(&mut [0.0; 2]).unwrap();
    assert!(DecimatedFrameRenderer::new(started, 6000, 2000, 5).is_err());
    let mut stream = DecimatedFrameRenderer::new(fresh(), 6000, 2000, 5).unwrap();
    assert!(stream.output_frames_for(101).is_err()); // no dropping last two frames
    assert!(stream.output_frames_for(102).is_err()); // no padding beyond history
    assert!(stream.validate_frame_count(u64::MAX).is_err());
    assert!(stream.validate_frame_rate(6000).is_err());
    assert_eq!(stream.output_frames_for(99).unwrap(), 33);
    for length in [0, 3, 12] {
        let mut output = vec![123.0; length];
        assert!(stream.frame_block(&mut output).is_err());
        assert_eq!(output, vec![123.0; length]);
        assert_eq!(stream.source().frames_rendered(), 0);
    }
    stream.frame_block(&mut [0.0; 10]).unwrap();
    assert_eq!(stream.frames_rendered(), 5);
    assert_eq!(stream.source().frames_rendered(), 15);
}

struct CancelOnPayload<'a> { cursor: Cursor<Vec<u8>>, gate: &'a CancelGate }
impl Write for CancelOnPayload<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let payload = self.cursor.position() >= 68;
        let written = self.cursor.write(bytes)?;
        if payload { self.gate.request(); }
        Ok(written)
    }
    fn flush(&mut self) -> std::io::Result<()> { self.cursor.flush() }
}
impl Seek for CancelOnPayload<'_> {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> { self.cursor.seek(position) }
}

#[test]
fn stream_cancellation_drains_whole_array_groups_and_keeps_filter_history_for_resume() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let make = || DecimatedFrameRenderer::new(
        array(&left, &right, &samples, 8000, 4*53, 1), 8000, 2000, 7,
    ).unwrap();
    let mut renderer = make();
    let gate = CancelGate::new_clock_free();
    let sink = CancelOnPayload { cursor: Cursor::new(Vec::new()), gate: &gate };
    let mut stream = Pcm16WavStream::new_interleaved(sink, 2000, 2, 0.02, 7).unwrap();
    assert_eq!(render_pressure_frames_pcm16(&mut renderer, &mut stream, &gate,
        &mut [0.0; 14], 53).unwrap(), WavFrameProgress::Cancelled { frames: 7 });
    assert_eq!(renderer.frames_rendered(), 7);
    assert_eq!(renderer.source().frames_rendered(), 28);
    assert_eq!(stream.frames_written(), 7);
    assert_eq!(render_pressure_frames_pcm16(&mut renderer, &mut stream,
        &CancelGate::new_clock_free(), &mut [0.0; 10], 46).unwrap(),
        WavFrameProgress::Completed { frames: 46 });
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(summary.frames, 53);
    let mut expected = make();
    let mut output = [0.0; 106];
    for part in output.chunks_mut(14) { expected.frame_block(part).unwrap(); }
    assert_eq!(sink.cursor.into_inner(), encode_pcm16_wav_interleaved(&output, 2000, 2, 0.02).unwrap().0);
}

// Fault injection only at the public producer boundary: the actual filter and
// adapter execute unchanged, including a source failure after internal progress.
struct FaultFrames { completed: u64, fail_at: u64, nonfinite: bool }
impl PressureFrameRenderer for FaultFrames {
    fn channel_count(&self) -> usize { 2 }
    fn frames_rendered(&self) -> u64 { self.completed }
    fn max_frame_block_len(&self) -> usize { 1 }
    fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate == 8000 { Ok(()) } else { Err(RenderError::Sizing { what: "fixture rate" }) }
    }
    fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        if frames <= 100-self.completed { Ok(()) } else { Err(RenderError::Sizing { what: "fixture horizon" }) }
    }
    fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        assert_eq!(output.len(), 2);
        let failed = self.completed == self.fail_at;
        self.completed += 1;
        if failed && !self.nonfinite { return Err(RenderError::Control { what: "late source failure" }); }
        output.copy_from_slice(&[0.001, if failed { f64::NAN } else { -0.002 }]);
        Ok(())
    }
}

#[test]
fn source_or_filter_failure_poison_the_array_instead_of_restarting_missing_history() {
    for nonfinite in [false, true] {
        let mut stream = DecimatedFrameRenderer::new(FaultFrames {
            completed: 0, fail_at: 4, nonfinite,
        }, 8000, 2000, 3).unwrap();
        assert!(stream.frame_block(&mut [0.0; 6]).is_err());
        let consumed = stream.source().frames_rendered();
        assert!(consumed >= 5);
        assert!(matches!(stream.frame_block(&mut [0.0; 2]), Err(RenderError::Poisoned)));
        assert!(matches!(stream.validate_frame_rate(2000), Err(RenderError::Poisoned)));
        assert!(matches!(stream.validate_frame_count(1), Err(RenderError::Poisoned)));
        assert_eq!(stream.source().frames_rendered(), consumed);
    }
    // Direct filter sanity: distinct input channels never share a delay line.
    let mut filter = Decimator::new(3, 2).unwrap();
    for _ in 0..200 {
        let frame = filter.preview(&[1.0, 0.0, 1.0, 0.0, 1.0, 0.0]).unwrap();
        assert_eq!(frame[1], 0.0);
        filter.commit();
    }
}
