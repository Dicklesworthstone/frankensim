//! G1/G3/G4: reflected stress fields on the original microphone/PCM clock.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, TransientQuadrupoleSource,
};
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::observation::HarmonicSampleClock;
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::quadrupole::{
    QuadrupoleFieldSample, QuadrupolePressureBudget, QuadrupolePressureConfig,
    QuadrupoleSourceWindow, TransientQuadrupoleRenderer,
};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::TAU;
use std::io::Cursor;

const FREQUENCY: f64 = 80.0;
const POSITIONS: [[f64; 3]; 2] = [[0.03, -0.02, 0.1], [-0.04, 0.01, 0.07]];
const RECEIVER: [f64; 3] = [0.4, 0.2, 0.35];
type Reader = fn(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 320.0).unwrap()
}
fn plane() -> IdealAcousticPlane {
    IdealAcousticPlane::new(
        [0.0; 3],
        [0.0, 0.0, 1.0],
        IdealPlaneBoundary::PressureRelease,
    )
    .unwrap()
}
fn observer() -> PreparedTransientQuadrupoleObserver {
    PreparedTransientQuadrupoleObserver::new(medium(), RECEIVER, &POSITIONS, 2)
        .unwrap()
        .with_orthogonal_planes(&[plane()], 4)
        .unwrap()
}
fn microphone() -> FirstOrder {
    FirstOrder::new(0.35, [0.6, 0.8, 0.0]).unwrap()
}
fn phasor(source: usize) -> [C64; 6] {
    let scale = if source == 0 { 0.01 } else { -0.006 };
    [
        C64::new(0.7, 0.2),
        C64::new(-0.3, 0.4),
        C64::new(0.9, -0.1),
        C64::new(0.23, -0.3),
        C64::new(-0.11, 0.3),
        C64::new(0.4, 0.2),
    ]
    .map(|q| q.scale(scale))
}
fn real_at(q: C64, time: f64) -> f64 {
    let (s, c) = (TAU * FREQUENCY * time).sin_cos();
    q.re * c + q.im * s
}
fn source_at(index: usize, time: Dd) -> Result<TransientQuadrupoleSource, AeroacError> {
    let q = phasor(index);
    let w = TAU * FREQUENCY;
    let time = time.to_f64();
    Ok(TransientQuadrupoleSource {
        stress_integral_n_m_s: q.map(|v| real_at(C64::new(-v.im / w, v.re / w), time)),
        stress_moment_n_m: q.map(|v| real_at(v, time)),
        stress_rate_n_m_s: q.map(|v| real_at(C64::new(w * v.im, -w * v.re), time)),
        stress_acceleration_n_m_s2: q.map(|v| real_at(v.scale(-w * w), time)),
    })
}
fn windows() -> [QuadrupoleSourceWindow; 2] {
    [QuadrupoleSourceWindow {
        source_time_interval_s: [Dd::from_f64(-0.1), Dd::from_f64(0.2)],
        source_bandwidth_hz: FREQUENCY,
    }; 2]
}
fn config() -> QuadrupolePressureConfig {
    QuadrupolePressureConfig {
        clock: HarmonicSampleClock::new(Dd::ZERO, 8_000).unwrap(),
        samples: 128,
        max_block: 128,
        receiver_bandwidth_hz: FREQUENCY,
    }
}
fn budget() -> QuadrupolePressureBudget {
    QuadrupolePressureBudget {
        maximum_sources: 2,
        maximum_output_samples: 128,
        maximum_source_path_samples: 512,
        maximum_scratch_bytes: 1 << 20,
    }
}
fn build<'a>(
    observer: &'a PreparedTransientQuadrupoleObserver,
    windows: &'a [QuadrupoleSourceWindow],
) -> TransientQuadrupoleRenderer<'a, Reader> {
    TransientQuadrupoleRenderer::new(
        observer,
        windows,
        source_at as Reader,
        microphone(),
        config(),
        budget(),
        &mut || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * scale.max(1e-20),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[test]
fn g1_reflected_transient_pcm_matches_original_harmonic_tensor_and_vector_fields() {
    let observer = observer();
    let windows = windows();
    let mut renderer = build(&observer, &windows);
    let harmonic = PreparedQuadrupoleObserver::new(medium(), FREQUENCY, RECEIVER, &POSITIONS, 2)
        .unwrap()
        .with_orthogonal_planes(&[plane()], 4)
        .unwrap()
        .evaluate(&[0, 1].map(|i| QuadrupolePhasor {
            stress_moment_n_m: phasor(i),
        }))
        .unwrap();
    let mut fields = [QuadrupoleFieldSample::default(); 128];
    renderer.field_block(&mut fields).unwrap();
    let velocity_scale = harmonic
        .particle_velocity_m_s
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max);
    for (index, field) in fields.iter().enumerate() {
        let time = renderer.receiver_time_at(index as u64).unwrap().to_f64();
        assert_eq!(field.receiver_time_s, time);
        close(
            field.pressure_pa,
            real_at(harmonic.pressure_pa, time),
            harmonic.pressure_pa.abs(),
        );
        for axis in 0..3 {
            close(
                field.particle_velocity_m_s[axis],
                real_at(harmonic.particle_velocity_m_s[axis], time),
                velocity_scale,
            );
            assert_eq!(
                field.instantaneous_intensity_w_m2[axis],
                field.pressure_pa * field.particle_velocity_m_s[axis]
            );
        }
        let directional = 0.35 * field.pressure_pa
            - 0.65
                * 1.2
                * 320.0
                * (0.6 * field.particle_velocity_m_s[0] + 0.8 * field.particle_velocity_m_s[1]);
        close(
            field.microphone_pressure_pa,
            directional,
            harmonic.pressure_pa.abs() + 1.2 * 320.0 * velocity_scale,
        );
    }
    assert!(renderer.block(&mut [0.0]).is_err());
}

#[test]
fn g3_fields_pcm_decimation_and_streaming_preserve_one_partition_invariant_clock() {
    let observer = observer();
    let windows = windows();
    let mut whole = build(&observer, &windows);
    let mut fields = [QuadrupoleFieldSample::default(); 128];
    whole.field_block(&mut fields).unwrap();
    let expected: Vec<_> = fields.iter().map(|f| f.microphone_pressure_pa).collect();
    let mut split = build(&observer, &windows);
    let mut first = [QuadrupoleFieldSample::default(); 7];
    split.field_block(&mut first).unwrap();
    assert_eq!(&first, &fields[..7]);
    let mut rest = vec![0.0; 121];
    for chunk in rest.chunks_mut(13) {
        PressureRenderer::block(&mut split, chunk).unwrap();
    }
    assert_eq!(rest, expected[7..]);
    assert_eq!(split.samples_rendered(), 128);

    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 8_000, 1.0, 17).unwrap();
    let result = render_pressure_pcm16(
        &mut build(&observer, &windows),
        &mut stream,
        &CancelGate::new_clock_free(),
        &mut [0.0; 17],
        128,
    )
    .unwrap();
    assert_eq!(result, ScheduledWavProgress::Completed { samples: 128 });
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(summary.samples, 128);
    assert_eq!(
        sink.into_inner(),
        encode_pcm16_wav(&expected, 8_000, 1.0).unwrap().0
    );

    let mut decimated =
        DecimatedRenderer::new(build(&observer, &windows), 8_000, 2_000, 32).unwrap();
    let mut one_block = [0.0; 32];
    decimated.block(&mut one_block).unwrap();
    let mut decimated_split =
        DecimatedRenderer::new(build(&observer, &windows), 8_000, 2_000, 32).unwrap();
    let mut several = [0.0; 32];
    for chunk in several.chunks_mut(7) {
        decimated_split.block(chunk).unwrap();
    }
    assert_eq!(one_block, several);
    assert!(decimated.info().delay_output_samples > 0.0);
}

#[test]
fn g3_large_epoch_translation_preserves_reflected_source_phase_on_the_original_dd_clock() {
    let observer = observer();
    let windows = windows();
    let mut expected = [0.0; 128];
    build(&observer, &windows).block(&mut expected).unwrap();
    let epoch = Dd::from_f64(2.0_f64.powi(40));
    let shifted_windows = windows.map(|window| QuadrupoleSourceWindow {
        source_time_interval_s: window.source_time_interval_s.map(|time| epoch + time),
        ..window
    });
    let shifted_config = QuadrupolePressureConfig {
        clock: HarmonicSampleClock::new(epoch, 8_000).unwrap(),
        ..config()
    };
    let mut shifted = TransientQuadrupoleRenderer::new(
        &observer,
        &shifted_windows,
        |index, time| source_at(index, time - epoch),
        microphone(),
        shifted_config,
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut actual = [0.0; 128];
    shifted.block(&mut actual).unwrap();
    let scale = expected.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    for (a, e) in actual.into_iter().zip(expected) {
        close(a, e, scale);
    }
    assert!(
        shifted
            .receiver_time_at(0)
            .unwrap()
            .lt(shifted.receiver_time_at(1).unwrap())
    );
}

#[test]
fn g4_refused_or_cancelled_callback_preserves_output_and_clock() {
    let observer = observer();
    let windows = windows();
    let mut renderer = build(&observer, &windows);
    let mut output = [123.0; 128];
    let mut polls = 0;
    assert!(
        renderer
            .block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 12
            })
            .is_err()
    );
    assert_eq!(output, [123.0; 128]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.block(&mut output).unwrap();
    let mut expected = [0.0; 128];
    build(&observer, &windows).block(&mut expected).unwrap();
    assert_eq!(output, expected);

    let refuse = |i, time| {
        if Dd::from_f64(0.005).lt(time) {
            Err(AeroacError::NonFinite {
                what: "test history refusal",
            })
        } else {
            source_at(i, time)
        }
    };
    let mut bad = TransientQuadrupoleRenderer::new(
        &observer,
        &windows,
        refuse,
        FirstOrder::default(),
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut untouched = [QuadrupoleFieldSample {
        pressure_pa: 17.0,
        ..Default::default()
    }; 128];
    let original = untouched;
    assert!(bad.field_block(&mut untouched).is_err());
    assert_eq!(untouched, original);
    assert_eq!(bad.samples_rendered(), 0);
}

#[test]
fn g0_complete_reflected_history_and_all_image_work_are_admitted_before_reading_sources() {
    let observer = observer();
    let mut windows = windows();
    let admit = |w: &[QuadrupoleSourceWindow], c, b| {
        TransientQuadrupoleRenderer::new(
            &observer,
            w,
            |_, _| -> Result<TransientQuadrupoleSource, AeroacError> {
                panic!("admission must not query history")
            },
            FirstOrder::default(),
            c,
            b,
            &mut || false,
        )
        .is_ok()
    };
    assert!(admit(&windows, config(), budget()));
    let [shortest, longest] = observer.source_delay_bounds_s(0).unwrap();
    assert!(longest > shortest);
    windows[0].source_time_interval_s[0] = Dd::from_f64(-0.5 * (shortest + longest));
    assert!(!admit(&windows, config(), budget()));
    windows[0].source_time_interval_s[0] = Dd::from_f64(-0.1);
    let mut work = budget();
    work.maximum_source_path_samples -= 1;
    assert!(!admit(&windows, config(), work));
    let mut memory = budget();
    memory.maximum_scratch_bytes = config().max_block * size_of::<QuadrupoleFieldSample>() - 1;
    assert!(!admit(&windows, config(), memory));
    let mut clock = config();
    clock.receiver_bandwidth_hz = 0.5 * f64::from(clock.clock.sample_rate_hz());
    assert!(!admit(&windows, clock, budget()));
    let mut bandwidth = config();
    bandwidth.receiver_bandwidth_hz = FREQUENCY - 1.0;
    assert!(!admit(&windows, bandwidth, budget()));
    let mut malformed = windows;
    malformed[0].source_time_interval_s[0] = Dd::from_pair(f64::MAX, f64::MAX);
    assert!(!admit(&malformed, config(), budget()));
}
