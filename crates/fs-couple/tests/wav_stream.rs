use std::io::{Cursor, Seek, SeekFrom, Write};
use std::sync::Arc;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::stream::{
    MAX_PCM16_WAV_SAMPLES, Pcm16WavStream, ScheduledWavProgress, WavStreamError,
    render_scheduled_pcm16,
};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn image() -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * 220.0,
        damping_ratio: 0.03,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn performance() -> ScheduledRenderer {
    let voice = ModalForceVoice::new(image(), vec![vec![2.0]], vec![1.0],
        ForceInitialization::RetainState).unwrap();
    ScheduledRenderer::from_modal_forces(vec![voice], vec![
        ModalForceEvent { sample: 37, voice: 0, port: 0, force_n: 0.0 },
        ModalForceEvent { sample: 71, voice: 0, port: 0, force_n: -0.5 },
    ], ForceRenderConfig { sample_rate_hz: 48_000, max_block: 128,
        max_events: 2, max_controls: 2, max_projection_terms: 3 }).unwrap()
}
fn reference() -> Vec<f64> {
    let mut model = image();
    (0..257).map(|n| model.step(&[if n < 37 { 2.0 } else if n < 71 { 0.0 } else { -1.0 }])
        .unwrap().observer_pressure_pa).collect()
}

#[test]
fn streaming_pcm_is_byte_exact_to_the_existing_encoder_across_partitions() {
    let mut pressure: Vec<_> = (0..257).map(|i| (i as f64-128.0)/73.0).collect();
    pressure.extend_from_slice(&[0.0, -0.0, 1.0, -1.0, f64::MAX, -f64::MAX]);
    let (expected, clips) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    for block in [1, 7, 37, 128, 1024] {
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 1.0, block).unwrap();
        for chunk in pressure.chunks(block) { stream.write_block(chunk).unwrap(); }
        let (cursor, report) = stream.finish().unwrap();
        assert_eq!(cursor.position(), expected.len() as u64);
        assert_eq!(cursor.into_inner(), expected);
        assert_eq!(report.samples, pressure.len() as u64);
        assert_eq!(report.clipped_samples, clips as u64);
        assert_eq!(report.full_scale_pa, 1.0);
    }
}

#[test]
fn invalid_blocks_do_not_change_payload_counters_or_poison_a_healthy_stream() {
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 1.0, 4).unwrap();
    for bad in [vec![0.1, f64::NAN], vec![0.1, f64::INFINITY], vec![0.0; 5]] {
        assert!(stream.write_block(&bad).is_err());
        assert_eq!(stream.samples_written(), 0);
        assert_eq!(stream.sink().get_ref().len(), 44);
    }
    stream.write_block(&[]).unwrap();
    stream.write_block(&[0.25, -0.25]).unwrap();
    let (cursor, report) = stream.finish().unwrap();
    assert_eq!(report.samples, 2);
    assert_eq!(cursor.into_inner(), encode_pcm16_wav(&[0.25,-0.25],48_000,1.0).unwrap().0);
}

#[test]
fn append_position_preserves_existing_bytes_and_interior_positions_refuse() {
    let mut existing = Cursor::new(b"existing".to_vec());
    assert!(Pcm16WavStream::new(&mut existing, 48_000, 1.0, 32).is_err());
    assert_eq!(existing.position(), 0);
    assert_eq!(existing.get_ref().as_slice(), b"existing");
    existing.set_position(8);
    let mut stream = Pcm16WavStream::new(existing, 48_000, 1.0, 32).unwrap();
    stream.write_block(&[0.25]).unwrap();
    let (cursor, _) = stream.finish().unwrap();
    let bytes = cursor.into_inner();
    assert_eq!(&bytes[..8], b"existing");
    assert_eq!(&bytes[8..], encode_pcm16_wav(&[0.25],48_000,1.0).unwrap().0.as_slice());
}

struct PartialFailure { cursor: Cursor<Vec<u8>>, remaining: usize }
impl Write for PartialFailure {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self.remaining == 0 { return Err(std::io::Error::other("injected partial write")); }
        let count = bytes.len().min(self.remaining);
        let n = self.cursor.write(&bytes[..count])?;
        self.remaining -= n;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> { Ok(()) }
}
impl Seek for PartialFailure {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> { self.cursor.seek(position) }
}

#[test]
fn partial_io_never_finalizes_or_resumes_as_a_complete_wav() {
    let output = PartialFailure { cursor: Cursor::new(Vec::new()), remaining: 50 };
    let mut stream = Pcm16WavStream::new(output,48_000,1.0,32).unwrap();
    assert!(matches!(stream.write_block(&[0.1; 8]), Err(WavStreamError::Io(_))));
    assert_eq!(stream.samples_written(), 0);
    assert!(matches!(stream.write_block(&[0.0]), Err(WavStreamError::Poisoned)));
    assert!(matches!(stream.finish(), Err(WavStreamError::Poisoned)));
}

#[test]
fn scheduled_force_export_equals_direct_physics_and_retains_the_last_short_block() {
    let (expected, _) = encode_pcm16_wav(&reference(),48_000,0.01).unwrap();
    for block in [1, 7, 37, 64, 128] {
        let mut renderer = performance();
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()),48_000,0.01,block).unwrap();
        let result = render_scheduled_pcm16(&mut renderer, &mut stream, &CancelGate::new(),
            &mut vec![0.0; block], 257).unwrap();
        assert_eq!(result, ScheduledWavProgress::Completed { samples: 257 });
        assert_eq!(renderer.samples_rendered(), 257);
        let (cursor, report) = stream.finish().unwrap();
        assert_eq!(report.samples, 257);
        assert_eq!(cursor.into_inner(), expected);
    }
}

struct CancelOnPayload { cursor: Cursor<Vec<u8>>, gate: Arc<CancelGate>, armed: bool }
impl Write for CancelOnPayload {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let n = self.cursor.write(bytes)?;
        if self.armed && self.cursor.position() > 44 { self.gate.request(); self.armed = false; }
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> { Ok(()) }
}
impl Seek for CancelOnPayload {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> { self.cursor.seek(position) }
}

#[test]
fn cancellation_drains_one_callback_then_resume_preserves_pending_release_and_file_bytes() {
    let gate = Arc::new(CancelGate::new());
    let output = CancelOnPayload { cursor: Cursor::new(Vec::new()), gate: gate.clone(), armed: true };
    let mut stream = Pcm16WavStream::new(output,48_000,0.01,37).unwrap();
    let mut renderer = performance();
    let mut scratch = [0.0; 37];
    assert_eq!(render_scheduled_pcm16(&mut renderer,&mut stream,&gate,&mut scratch,257).unwrap(),
        ScheduledWavProgress::Cancelled { samples: 37 });
    assert_eq!(renderer.samples_rendered(),37);
    assert_eq!(stream.samples_written(),37);
    assert_eq!(renderer.pending_controls()[0].sample,37);
    assert_eq!(render_scheduled_pcm16(&mut renderer,&mut stream,&CancelGate::new(),&mut scratch,220).unwrap(),
        ScheduledWavProgress::Completed { samples: 220 });
    let (output, report) = stream.finish().unwrap();
    assert_eq!(report.samples,257);
    assert_eq!(output.cursor.into_inner(), encode_pcm16_wav(&reference(),48_000,0.01).unwrap().0);
}

#[test]
fn rate_timeline_capacity_and_whole_request_overflow_refuse_before_physics() {
    let mut renderer = performance();
    let mut wrong_rate = Pcm16WavStream::new(Cursor::new(Vec::new()),44_100,1.0,128).unwrap();
    assert!(render_scheduled_pcm16(&mut renderer,&mut wrong_rate,&CancelGate::new(),&mut [0.0; 32],32).is_err());
    assert_eq!(renderer.samples_rendered(),0);
    assert_eq!(wrong_rate.samples_written(),0);
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()),48_000,1.0,128).unwrap();
    assert!(render_scheduled_pcm16(&mut renderer,&mut stream,&CancelGate::new(),&mut [0.0; 129],32).is_err());
    assert!(render_scheduled_pcm16(&mut renderer,&mut stream,&CancelGate::new(),&mut [0.0; 32],MAX_PCM16_WAV_SAMPLES+1).is_err());
    assert_eq!(renderer.samples_rendered(),0);
    assert_eq!(stream.samples_written(),0);
    stream.write_block(&[0.0]).unwrap();
    assert!(render_scheduled_pcm16(&mut renderer,&mut stream,&CancelGate::new(),&mut [0.0; 32],32).is_err());
    assert_eq!(renderer.samples_rendered(),0);
    assert_eq!(stream.samples_written(),1);
}

#[test]
fn already_cancelled_export_can_finalize_an_explicitly_empty_wav() {
    let gate = CancelGate::new(); gate.request();
    let mut renderer = performance();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()),48_000,1.0,32).unwrap();
    assert_eq!(render_scheduled_pcm16(&mut renderer,&mut stream,&gate,&mut [0.0;32],257).unwrap(),
        ScheduledWavProgress::Cancelled { samples: 0 });
    let (cursor, report) = stream.finish().unwrap();
    assert_eq!(report.samples,0);
    let bytes = cursor.into_inner();
    assert_eq!(bytes.len(),44);
    assert_eq!(&bytes[4..8],&36_u32.to_le_bytes());
    assert_eq!(&bytes[40..44],&0_u32.to_le_bytes());
}
