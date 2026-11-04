//! G1/G3: physical harmonic receivers use the original pressure/WAV owner.

use fs_aeroac::observation::{
    HarmonicFieldLine, HarmonicObservationBudget, HarmonicSampleBudget, HarmonicSampleClock,
    PreparedHarmonicField,
};
use fs_couple::pcm_wav::{
    encode_pcm16_wav,
    observation::{PressureRenderer, harmonic::HarmonicPressureRenderer},
    stream::{Pcm16WavStream, render_pressure_pcm16},
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, dd::Dd};
use std::io::Cursor;

fn field() -> PreparedHarmonicField<3> {
    PreparedHarmonicField::new(
        &[
            HarmonicFieldLine {
                frequency_hz: 213.17,
                pressure_pa: C64::new(0.2, -0.1),
                particle_velocity_m_s: [C64::ZERO; 3],
            },
            HarmonicFieldLine {
                frequency_hz: 217.83,
                pressure_pa: C64::new(0.07, 0.03),
                particle_velocity_m_s: [C64::ZERO; 3],
            },
        ],
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap()
}

fn budget() -> HarmonicSampleBudget {
    HarmonicSampleBudget {
        max_samples: 128,
        max_tone_samples: 256,
    }
}

#[test]
fn g1_harmonic_pressure_reaches_original_pcm_stream_without_gain_change() {
    let field = field();
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.031), 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(&field, clock, budget()).unwrap();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 8000, 1.0, 128).unwrap();
    let gate = CancelGate::new_clock_free();
    render_pressure_pcm16(&mut renderer, &mut stream, &gate, &mut [0.0; 37], 513).unwrap();
    let (output, summary) = stream.finish().unwrap();
    let mut pressure = vec![0.0; 513];
    field
        .sample_pressure_into(
            clock,
            0,
            &mut pressure,
            HarmonicSampleBudget {
                max_samples: 513,
                max_tone_samples: 1026,
            },
        )
        .unwrap();
    let (expected, clips) = encode_pcm16_wav(&pressure, 8000, 1.0).unwrap();
    assert_eq!(output.into_inner(), expected);
    assert_eq!(summary.samples, 513);
    assert_eq!(summary.clipped_samples, clips as u64);
    assert_eq!(renderer.samples_rendered(), 513);
}

#[test]
fn g3_harmonic_callback_partition_and_admission_are_transactional() {
    let field = field();
    let clock = HarmonicSampleClock::new(Dd::ZERO, 8000).unwrap();
    let mut a = HarmonicPressureRenderer::new(&field, clock, budget()).unwrap();
    let mut b = HarmonicPressureRenderer::new(&field, clock, budget()).unwrap();
    let mut whole = [0.0; 128];
    let mut pieces = [0.0; 128];
    a.block(&mut whole).unwrap();
    for chunk in pieces.chunks_mut(13) {
        b.block(chunk).unwrap();
    }
    assert_eq!(whole, pieces);
    let before = b.samples_rendered();
    let mut oversized = [17.0; 129];
    assert!(b.block(&mut oversized).is_err());
    assert_eq!(oversized, [17.0; 129]);
    assert_eq!(b.samples_rendered(), before);
    assert!(b.validate_sample_rate(44_100).is_err());
    let aliased = HarmonicSampleClock::new(Dd::ZERO, 400).unwrap();
    assert!(HarmonicPressureRenderer::new(&field, aliased, budget()).is_err());
}

#[test]
fn g0_unrepresentable_pressure_does_not_publish_a_partial_callback() {
    let field = PreparedHarmonicField::<1>::new(
        &[
            HarmonicFieldLine {
                frequency_hz: 100.0,
                pressure_pa: C64::new(1e308, 0.0),
                particle_velocity_m_s: [C64::ZERO],
            },
            HarmonicFieldLine {
                frequency_hz: 200.0,
                pressure_pa: C64::new(1e308, 0.0),
                particle_velocity_m_s: [C64::ZERO],
            },
        ],
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap();
    let clock = HarmonicSampleClock::new(Dd::ZERO, 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(&field, clock, budget()).unwrap();
    let mut output = [17.0; 8];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [17.0; 8]);
    assert_eq!(renderer.samples_rendered(), 0);
}
