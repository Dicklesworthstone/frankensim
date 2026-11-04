//! Simultaneous physical interference reaches the original calibrated WAV owner.
use fs_aeroac::observation::{
    HarmonicEnsembleBudget, HarmonicFactorStatistics, HarmonicFieldFactor, HarmonicFieldLine,
    HarmonicSampleBudget, HarmonicSampleClock, PreparedHarmonicEnsemble,
    array::{
        HarmonicArrayBudget, HarmonicArrayReceiver, HarmonicFactorCoordinate, HarmonicFactorDraw,
        HarmonicReadoutBudget, HarmonicReadoutTap, HarmonicRealizationBudget,
        PreparedHarmonicArray, PreparedHarmonicReadout,
    },
};
use fs_couple::pcm_wav::{
    encode_pcm16_wav_interleaved,
    observation::{
        PressureFrameRenderer,
        harmonic::{HarmonicOutputChannel, HarmonicOutputSource, HarmonicPressureFrameRenderer},
    },
    stream::{Pcm16WavStream, WavFrameProgress, WavStreamError, render_pressure_frames_pcm16},
};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::io::{Cursor, Seek, SeekFrom, Write};

const COORDINATES: [HarmonicFactorCoordinate<'static>; 2] = [
    HarmonicFactorCoordinate {
        id: "common real",
        statistics: HarmonicFactorStatistics::Real,
    },
    HarmonicFactorCoordinate {
        id: "common proper",
        statistics: HarmonicFactorStatistics::ProperComplex,
    },
];
const PREPARE: HarmonicRealizationBudget = HarmonicRealizationBudget {
    max_channels: 4,
    max_factors: 2,
    max_total_lines: 64,
    max_line_pairs: 4096,
};
const SAMPLE: HarmonicSampleBudget = HarmonicSampleBudget {
    max_samples: 512,
    max_tone_samples: 2048,
};
const READOUT: HarmonicReadoutBudget = HarmonicReadoutBudget {
    max_taps: 2,
    max_total_lines: 32,
    max_line_pairs: 1024,
};
const DRAWS: [HarmonicFactorDraw; 2] = [
    HarmonicFactorDraw::Real(-0.5),
    HarmonicFactorDraw::ProperComplex(C64::new(0.4, -0.6)),
];

fn bank(scale: f64) -> PreparedHarmonicEnsemble<3> {
    let line = |f, p, u| HarmonicFieldLine {
        frequency_hz: f,
        pressure_pa: C64::new(p * scale, 0.0),
        particle_velocity_m_s: [C64::new(0.0, u * scale), C64::ZERO, C64::ZERO],
    };
    PreparedHarmonicEnsemble::new(
        &[line(80.0, 0.3, 0.021), line(130.0, -0.11, 0.013)],
        &[
            HarmonicFieldFactor {
                lines: &[line(100.0, 0.25, -0.04)],
                statistics: HarmonicFactorStatistics::Real,
            },
            HarmonicFieldFactor {
                lines: &[line(160.0, 0.17, 0.035)],
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ],
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: 4,
            max_line_pairs: 6,
        },
    )
    .unwrap()
}
fn array<'a>(
    left: &'a PreparedHarmonicEnsemble<3>,
    right: &'a PreparedHarmonicEnsemble<3>,
) -> PreparedHarmonicArray<'a, 3> {
    let medium = FluidAcousticState::try_new(2.0, 3.0).unwrap();
    PreparedHarmonicArray::new(
        "physical source clock",
        &COORDINATES,
        &[
            HarmonicArrayReceiver {
                position_m: [0.0; 3],
                acoustic: medium,
                ensemble: left,
            },
            HarmonicArrayReceiver {
                position_m: [0.3, 0.0, 0.0],
                acoustic: medium,
                ensemble: right,
            },
        ],
        HarmonicArrayBudget {
            max_receivers: 2,
            max_factors: 2,
            max_total_lines: 8,
            max_line_pairs: 64,
        },
    )
    .unwrap()
}
fn readouts<'a>(
    array: &PreparedHarmonicArray<'a, 3>,
) -> (PreparedHarmonicReadout<'a>, PreparedHarmonicReadout<'a>) {
    let directional = array
        .prepare_readout(
            &[HarmonicReadoutTap {
                receiver: 0,
                gain: -0.5,
                delay_s: 0.0013,
                pressure_fraction: 0.25,
                front_axis: [1.0, 0.0, 0.0],
            }],
            READOUT,
        )
        .unwrap();
    let null = array
        .prepare_readout(
            &[
                HarmonicReadoutTap::pressure(0, 2.0, 0.0),
                HarmonicReadoutTap::pressure(1, 1.0, 0.0),
            ],
            READOUT,
        )
        .unwrap();
    (directional, null)
}
fn channels<'a>(
    directional: &'a PreparedHarmonicReadout<'a>,
    null: &'a PreparedHarmonicReadout<'a>,
) -> [HarmonicOutputChannel<'a>; 4] {
    [
        HarmonicOutputChannel {
            id: "receiver right",
            source: HarmonicOutputSource::Receiver(1),
        },
        HarmonicOutputChannel {
            id: "directional delayed",
            source: HarmonicOutputSource::Readout(directional),
        },
        HarmonicOutputChannel {
            id: "receiver left",
            source: HarmonicOutputSource::Receiver(0),
        },
        HarmonicOutputChannel {
            id: "coherent null",
            source: HarmonicOutputSource::Readout(null),
        },
    ]
}
fn clock() -> HarmonicSampleClock {
    HarmonicSampleClock::new(Dd::from_f64(0.031), 8_000).unwrap()
}

#[test]
fn original_shared_draw_directional_velocity_and_delay_reach_simultaneous_channels() {
    let (left, right) = (bank(1.0), bank(-2.0));
    let array = array(&left, &right);
    let realized = array.realize(&DRAWS, PREPARE).unwrap();
    let (directional, null) = readouts(&array);
    let channels = channels(&directional, &null);
    let mut whole =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let mut split =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let mut p = [0.0; 512];
    let mut chunks = [0.0; 512];
    whole.frame_block(&mut p).unwrap();
    for block in chunks.chunks_mut(4 * 13) {
        split.frame_block(block).unwrap();
    }
    assert_eq!(p, chunks);
    assert_eq!(whole.frames_rendered(), 128);
    assert_eq!(whole.channels()[1].id, "directional delayed");
    assert_eq!(whole.realization().factor_draws(), &DRAWS);
    for (i, frame) in p.chunks_exact(4).enumerate() {
        let time = clock().time_at(i as u64).unwrap();
        let physical = realized.receivers()[0]
            .field
            .field_at_clock_dd(time)
            .unwrap();
        let delayed = realized.receivers()[0]
            .field
            .field_at_clock_dd(time - Dd::from_f64(0.0013))
            .unwrap();
        let expected =
            -0.5 * (0.25 * delayed.pressure_pa - 0.75 * 6.0 * delayed.particle_velocity_m_s[0]);
        assert_eq!(frame[2], physical.pressure_pa);
        assert_eq!(frame[0], -2.0 * frame[2]);
        assert!((frame[1] - expected).abs() < 3e-15);
        assert_eq!(frame[3], 0.0);
    }
    let (expected, clips) = encode_pcm16_wav_interleaved(&p, 8_000, 4, 2.0).unwrap();
    for block_frames in [1, 7, 37, 128] {
        let mut renderer =
            HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE)
                .unwrap();
        let mut stream =
            Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 8_000, 4, 2.0, 128).unwrap();
        assert_eq!(
            render_pressure_frames_pcm16(
                &mut renderer,
                &mut stream,
                &CancelGate::new_clock_free(),
                &mut vec![0.0; block_frames * 4],
                128
            )
            .unwrap(),
            WavFrameProgress::Completed { frames: 128 }
        );
        let (out, summary) = stream.finish().unwrap();
        assert_eq!(out.into_inner(), expected);
        assert_eq!(
            (
                summary.channels,
                summary.frames,
                summary.samples,
                summary.clipped_samples
            ),
            (4, 128, 512, clips as u64)
        );
    }
}

#[test]
fn shape_budgets_dark_nyquist_and_late_cancellation_publish_no_partial_frames() {
    let (left, right) = (bank(1.0), bank(-2.0));
    let array = array(&left, &right);
    let realized = array.realize(&DRAWS, PREPARE).unwrap();
    let (directional, null) = readouts(&array);
    let channels = channels(&directional, &null);
    let mut renderer =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let mut bad = [17.0; 7];
    assert!(renderer.frame_block(&mut bad).is_err());
    assert_eq!(bad, [17.0; 7]);
    let mut output = [17.0; 128];
    let mut polls = 0;
    assert!(
        renderer
            .frame_block_cancellable(&mut output, || {
                polls += 1;
                polls == 40
            })
            .is_err()
    );
    assert_eq!(polls, 40);
    assert_eq!(output, [17.0; 128]);
    assert_eq!(renderer.frames_rendered(), 0);
    renderer.frame_block(&mut output).unwrap();
    assert_eq!(renderer.frames_rendered(), 32);
    let mut short_budget = PREPARE;
    short_budget.max_total_lines = 12; // the four selected banks require sixteen
    assert!(
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), short_budget, SAMPLE)
            .is_err()
    );
    let mut sample_budget = SAMPLE;
    sample_budget.max_samples = 511;
    assert!(
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, sample_budget)
            .is_err()
    );
    sample_budget = SAMPLE;
    sample_budget.max_tone_samples = 2047;
    assert!(
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, sample_budget)
            .is_err()
    );
    let mut duplicate = channels;
    duplicate[1].id = duplicate[0].id;
    assert!(
        HarmonicPressureFrameRenderer::new(&realized, &duplicate, clock(), PREPARE, SAMPLE)
            .is_err()
    );
    // Exact physical cancellation does not erase the original out-of-band line.
    let aliased = HarmonicSampleClock::new(Dd::ZERO, 320).unwrap();
    assert!(
        HarmonicPressureFrameRenderer::new(&realized, &channels[3..], aliased, PREPARE, SAMPLE)
            .is_err()
    );
}

struct CancelOnPayload<'a> {
    cursor: Cursor<Vec<u8>>,
    gate: &'a CancelGate,
}
impl Write for CancelOnPayload<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let payload = self.cursor.position() >= 68;
        let n = self.cursor.write(bytes)?;
        if payload {
            self.gate.request();
        }
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.cursor.flush()
    }
}
impl Seek for CancelOnPayload<'_> {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
        self.cursor.seek(position)
    }
}

#[test]
fn cancellation_drains_all_channels_and_continuation_keeps_the_same_frame_clock() {
    let (left, right) = (bank(1.0), bank(-2.0));
    let array = array(&left, &right);
    let realized = array.realize(&DRAWS, PREPARE).unwrap();
    let (directional, null) = readouts(&array);
    let channels = channels(&directional, &null);
    let mut renderer =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let gate = CancelGate::new_clock_free();
    let sink = CancelOnPayload {
        cursor: Cursor::new(Vec::new()),
        gate: &gate,
    };
    let mut stream = Pcm16WavStream::new_interleaved(sink, 8_000, 4, 2.0, 16).unwrap();
    assert_eq!(
        render_pressure_frames_pcm16(&mut renderer, &mut stream, &gate, &mut [0.0; 28], 23)
            .unwrap(),
        WavFrameProgress::Cancelled { frames: 7 }
    );
    assert_eq!(
        (
            renderer.frames_rendered(),
            stream.frames_written(),
            stream.samples_written()
        ),
        (7, 7, 28)
    );
    assert_eq!(
        render_pressure_frames_pcm16(
            &mut renderer,
            &mut stream,
            &CancelGate::new_clock_free(),
            &mut [0.0; 16],
            16
        )
        .unwrap(),
        WavFrameProgress::Completed { frames: 16 }
    );
    let (sink, summary) = stream.finish().unwrap();
    let mut reference =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let mut p = [0.0; 92];
    reference.frame_block(&mut p).unwrap();
    assert_eq!(
        sink.cursor.into_inner(),
        encode_pcm16_wav_interleaved(&p, 8_000, 4, 2.0).unwrap().0
    );
    assert_eq!((summary.frames, summary.samples), (23, 92));
}

#[test]
fn stream_whole_request_admission_precedes_producer_or_sink_advancement() {
    let (left, right) = (bank(1.0), bank(-2.0));
    let array = array(&left, &right);
    let realized = array.realize(&DRAWS, PREPARE).unwrap();
    let (directional, null) = readouts(&array);
    let channels = channels(&directional, &null);
    let mut renderer =
        HarmonicPressureFrameRenderer::new(&realized, &channels, clock(), PREPARE, SAMPLE).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut stream =
        Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 8_000, 4, 2.0, 128).unwrap();
    let before = stream.sink().get_ref().clone();
    for count in [u64::MAX, u64::from(u32::MAX) / 8] {
        assert!(
            render_pressure_frames_pcm16(&mut renderer, &mut stream, &gate, &mut [0.0; 16], count)
                .is_err()
        );
    }
    assert!(
        render_pressure_frames_pcm16(&mut renderer, &mut stream, &gate, &mut [0.0; 15], 1).is_err()
    );
    assert!(stream.write_block(&[0.0; 3]).is_err());
    assert!(stream.write_block(&[0.0, 0.0, f64::NAN, 0.0]).is_err());
    assert_eq!(stream.sink().get_ref(), &before);
    assert_eq!(
        (renderer.frames_rendered(), stream.frames_written()),
        (0, 0)
    );
    let mut wrong_channels = Pcm16WavStream::new(Cursor::new(Vec::new()), 8_000, 2.0, 128).unwrap();
    assert!(
        render_pressure_frames_pcm16(&mut renderer, &mut wrong_channels, &gate, &mut [0.0; 16], 1)
            .is_err()
    );
    let mut wrong_rate =
        Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 16_000, 4, 2.0, 128).unwrap();
    assert!(
        render_pressure_frames_pcm16(&mut renderer, &mut wrong_rate, &gate, &mut [0.0; 16], 1)
            .is_err()
    );
    stream.write_block(&[0.0; 4]).unwrap();
    assert!(
        render_pressure_frames_pcm16(&mut renderer, &mut stream, &gate, &mut [0.0; 16], 1).is_err()
    );
    assert_eq!(renderer.frames_rendered(), 0);
}

struct PartialSink(Cursor<Vec<u8>>);
impl Write for PartialSink {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let remaining = 72usize.saturating_sub(self.0.position() as usize);
        if remaining == 0 {
            return Err(std::io::Error::other("partial payload"));
        }
        self.0.write(&bytes[..bytes.len().min(remaining)])
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
impl Seek for PartialSink {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
        self.0.seek(position)
    }
}

#[test]
fn incomplete_io_poisoning_never_reports_a_partial_physical_frame_as_complete() {
    let mut stream =
        Pcm16WavStream::new_interleaved(PartialSink(Cursor::new(Vec::new())), 8_000, 4, 2.0, 8)
            .unwrap();
    assert!(matches!(
        stream.write_block(&[0.25; 4]),
        Err(WavStreamError::Io(_))
    ));
    assert_eq!((stream.frames_written(), stream.samples_written()), (0, 0));
    assert_eq!(stream.sink().0.get_ref().len(), 72);
    assert!(matches!(
        stream.write_block(&[0.0; 4]),
        Err(WavStreamError::Poisoned)
    ));
    assert!(matches!(stream.finish(), Err(WavStreamError::Poisoned)));
}
