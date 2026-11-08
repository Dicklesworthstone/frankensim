//! Original incident-plus-scattered BEM harmonics reach the existing PCM owner.
use fs_bem::helmholtz::Medium;
use fs_bem::near_field::{FirstOrder, Options};
use fs_bem::observation_transfers::{
    HarmonicExposureFactorStatistics, HarmonicExposureLimits, HarmonicExposureLine,
    HarmonicFactorDraw, HarmonicField, HarmonicFieldSample, HarmonicPressureArrayLimits,
    HarmonicPressureChannel, HarmonicPressureSampleLimits, HarmonicPressureTap,
    HarmonicRealizationLimits, HarmonicSampleClock, HarmonicSampleLimits, PointSourceSpec,
    PreparedHarmonicExposure, PreparedHarmonicObservation, PreparedHarmonicPressureSignal,
    PreparedHarmonicRealization, ReceiverSpec, ScatteringObservation, TransferError,
};
use fs_bem::panel3d::SpherePanels;
use fs_couple::pcm_wav::{
    encode_pcm16_wav_interleaved,
    observation::{PressureFrameRenderer, harmonic::BemHarmonicPressureFrameRenderer},
    stream::{Pcm16WavStream, WavFrameProgress, render_pressure_frames_pcm16},
};
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_math::{c64::C64, dd::Dd};
use std::io::Cursor;

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 343.0,
};
const OMEGAS: [f64; 2] = [
    core::f64::consts::TAU * 80.0,
    core::f64::consts::TAU * 130.0,
];
const MEAN: [C64; 2] = [C64::new(2e-5, -1e-5), C64::new(-8e-6, 1.3e-5)];
const REAL: [C64; 2] = [C64::new(7e-6, 2e-6), C64::new(-2e-6, 4e-6)];
const PROPER: [C64; 2] = [C64::new(-3e-6, 5e-6), C64::new(4e-6, -2e-6)];
const SAMPLE: HarmonicPressureSampleLimits = HarmonicPressureSampleLimits {
    maximum_frames: 128,
    maximum_channel_line_samples: 1024,
};

fn owners() -> [PreparedHarmonicObservation; 2] {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let sources = [PointSourceSpec {
        id: "prescribed-volume-flow".into(),
        position_m: [0.08, 0.02, 0.03],
    }];
    let receivers = [
        ReceiverSpec {
            id: "near".into(),
            position_m: [0.04, 0.05, 0.02],
        },
        ReceiverSpec {
            id: "shadow".into(),
            position_m: [-0.06, 0.03, 0.04],
        },
    ];
    let observation = ScatteringObservation::new(
        &surface,
        MEDIUM,
        &sources,
        &receivers,
        core::f64::consts::TAU * 300.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap();
    OMEGAS.map(|omega| observation.prepare_harmonic(omega, || false).unwrap())
}

fn realization(owners: &[PreparedHarmonicObservation; 2]) -> PreparedHarmonicRealization<'_> {
    let lines = |a: &'static [C64; 2]| {
        [
            HarmonicExposureLine {
                observation: &owners[0],
                coefficients: &a[..1],
            },
            HarmonicExposureLine {
                observation: &owners[1],
                coefficients: &a[1..],
            },
        ]
    };
    PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "physical-volume-flow-clock",
        &lines(&MEAN),
        &[&lines(&REAL), &lines(&PROPER)],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        HarmonicExposureLimits {
            maximum_lines: 8,
            maximum_factors: 2,
            maximum_transfer_products: 1024,
            maximum_receiver_frequency_pairs: 1024,
            maximum_workspace_bytes: 1_000_000,
        },
    )
    .unwrap()
    .realize(
        &[
            HarmonicFactorDraw::Real(-0.7),
            HarmonicFactorDraw::ProperComplex(C64::new(0.3, -1.1)),
        ],
        HarmonicRealizationLimits {
            maximum_lines: 8,
            maximum_field_products: 1024,
            maximum_workspace_bytes: 1_000_000,
        },
    )
    .unwrap()
}

fn prepare() -> HarmonicPressureArrayLimits {
    HarmonicPressureArrayLimits {
        maximum_channels: 4,
        maximum_taps: 16,
        maximum_prepare_products: 8192,
        maximum_query_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn clock() -> HarmonicSampleClock {
    HarmonicSampleClock::new(Dd::from_f64(0.031), 8000).unwrap()
}
fn patterns() -> [FirstOrder; 2] {
    [
        FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.0, [0.0, 1.0, 0.0]).unwrap(),
    ]
}
fn tap(receiver: usize, weight: f64, delay_s: f64) -> HarmonicPressureTap {
    HarmonicPressureTap {
        receiver,
        weight,
        delay_s,
    }
}
fn channels() -> [HarmonicPressureChannel; 4] {
    [
        HarmonicPressureChannel {
            id: "near-cardioid".into(),
            taps: vec![tap(0, 1.0, 0.0)],
        },
        HarmonicPressureChannel {
            id: "signed-beam".into(),
            taps: vec![tap(0, 0.7, 0.00013), tap(1, -0.4, 0.00037)],
        },
        HarmonicPressureChannel {
            id: "shadow-dipole".into(),
            taps: vec![tap(1, 1.0, 0.00021)],
        },
        HarmonicPressureChannel {
            id: "coherent-null".into(),
            taps: vec![tap(0, 1.0, 0.0), tap(0, -1.0, 0.0)],
        },
    ]
}

#[test]
fn scattered_directional_harmonics_reach_original_multichannel_wav_with_the_same_draw() {
    let owners = owners();
    let realization = realization(&owners);
    let channels = channels();
    let patterns = patterns();
    let signal =
        PreparedHarmonicPressureSignal::new(&realization, &channels, Some(&patterns), prepare())
            .unwrap();
    let mut whole = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    let mut split = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    let mut pressure = [0.0; 512];
    let mut blocks = [0.0; 512];
    whole.frame_block(&mut pressure).unwrap();
    for block in blocks.chunks_mut(4 * 13) {
        split.frame_block(block).unwrap();
    }
    assert_eq!(pressure, blocks);
    assert_eq!(whole.channels(), &channels);
    assert_eq!(whole.frames_rendered(), 128);
    let mut fields = [[HarmonicField::default(); 2]; 2];
    for f in 0..2 {
        let coefficient = MEAN[f] + REAL[f].scale(-0.7) + PROPER[f] * C64::new(0.3, -1.1);
        owners[f]
            .evaluate_coherent(&[coefficient], &mut fields[f], || false)
            .unwrap();
    }
    for (frame, samples) in pressure.chunks_exact(4).enumerate() {
        let time = 0.031 + frame as f64 / 8000.0;
        for (channel, &sample) in channels.iter().zip(samples) {
            let mut expected = 0.0;
            for tap in &channel.taps {
                for f in 0..2 {
                    let field = fields[f][tap.receiver];
                    let angle = OMEGAS[f] * (time - tap.delay_s);
                    let (sin, cos) = angle.sin_cos();
                    let real = |a: C64| a.re * cos + a.im * sin;
                    let p = real(field.pressure_pa);
                    let velocity = field.particle_velocity_m_s.map(real);
                    let pattern = patterns[tap.receiver];
                    let axial: f64 = pattern
                        .front_axis()
                        .iter()
                        .zip(velocity)
                        .map(|(n, u)| n * u)
                        .sum();
                    expected += tap.weight
                        * (pattern.pressure_fraction() * p
                            - (1.0 - pattern.pressure_fraction())
                                * MEDIUM.density
                                * MEDIUM.sound_speed
                                * axial);
                }
            }
            assert!(
                (sample - expected).abs() < 5e-12,
                "{sample:e} != {expected:e}"
            );
        }
        assert_eq!(samples[3], 0.0);
    }
    assert!(pressure.iter().any(|v| v.abs() > 1e-4));
    let (expected, clips) = encode_pcm16_wav_interleaved(&pressure, 8000, 4, 10.0).unwrap();
    for block_frames in [1, 7, 37, 128] {
        let mut renderer = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
        let mut stream =
            Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 8000, 4, 10.0, 128).unwrap();
        assert_eq!(
            render_pressure_frames_pcm16(
                &mut renderer,
                &mut stream,
                &CancelGate::new_clock_free(),
                &mut vec![0.0; 4 * block_frames],
                128,
            )
            .unwrap(),
            WavFrameProgress::Completed { frames: 128 }
        );
        let (bytes, summary) = stream.finish().unwrap();
        assert_eq!(bytes.into_inner(), expected);
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
fn callback_refusal_and_late_cancellation_publish_no_partial_channels_and_allow_retry() {
    let owners = owners();
    let realization = realization(&owners);
    let signal = PreparedHarmonicPressureSignal::new(
        &realization,
        &channels(),
        Some(&patterns()),
        prepare(),
    )
    .unwrap();
    let mut renderer = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    let mut wrong_shape = [17.0; 3];
    assert!(renderer.frame_block(&mut wrong_shape).is_err());
    assert_eq!(wrong_shape, [17.0; 3]);
    let mut oversized = [17.0; 516];
    assert!(renderer.frame_block(&mut oversized).is_err());
    assert_eq!(oversized, [17.0; 516]);
    assert!(renderer.validate_frame_rate(48_000).is_err());
    let mut output = [17.0; 64];
    let mut polls = 0;
    assert!(matches!(
        renderer.frame_block_cancellable(&mut output, || {
            polls += 1;
            polls == 20
        }),
        Err(RenderError::BemObservation(TransferError::Cancelled))
    ));
    assert_eq!(output, [17.0; 64]);
    assert_eq!(renderer.frames_rendered(), 0);
    renderer.frame_block(&mut output).unwrap();
    let mut expected = [0.0; 64];
    let mut fresh = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    fresh.frame_block(&mut expected).unwrap();
    assert_eq!(output, expected);
    assert_eq!(renderer.frames_rendered(), 16);
    assert!(renderer.validate_frame_count(u64::MAX).is_err());
    assert_eq!(renderer.frames_rendered(), 16);

    let dark = [HarmonicPressureChannel {
        id: "dark".into(),
        taps: vec![],
    }];
    let null = PreparedHarmonicPressureSignal::new(&realization, &dark, None, prepare()).unwrap();
    assert!(
        BemHarmonicPressureFrameRenderer::new(
            &null,
            HarmonicSampleClock::new(Dd::ZERO, 100).unwrap(),
            SAMPLE
        )
        .is_err()
    );
    assert!(
        BemHarmonicPressureFrameRenderer::new(
            &signal,
            clock(),
            HarmonicPressureSampleLimits {
                maximum_channel_line_samples: 0,
                ..SAMPLE
            }
        )
        .is_err()
    );
    let empty = PreparedHarmonicPressureSignal::new(&realization, &[], None, prepare()).unwrap();
    assert!(BemHarmonicPressureFrameRenderer::new(&empty, clock(), SAMPLE).is_err());
}

#[test]
fn pressure_only_channels_preserve_original_field_samples_and_wrong_sink_clock_refuses_early() {
    let owners = owners();
    let realization = realization(&owners);
    let channels: Vec<_> = (0..2)
        .map(|i| HarmonicPressureChannel {
            id: format!("physical-{i}"),
            taps: vec![tap(i, 1.0, 0.0)],
        })
        .collect();
    let signal =
        PreparedHarmonicPressureSignal::new(&realization, &channels, None, prepare()).unwrap();
    let mut renderer = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    let mut output = [0.0; 64];
    renderer.frame_block(&mut output).unwrap();
    let mut fields = [HarmonicFieldSample::default(); 64];
    realization
        .sample_fields_into(
            clock(),
            0,
            &mut fields,
            HarmonicSampleLimits {
                maximum_frames: 32,
                maximum_receiver_line_samples: 128,
            },
        )
        .unwrap();
    for (pressure, field) in output.iter().zip(&fields) {
        assert_eq!(pressure.to_bits(), field.pressure_pa.to_bits());
    }
    let mut renderer = BemHarmonicPressureFrameRenderer::new(&signal, clock(), SAMPLE).unwrap();
    let mut stream =
        Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 48_000, 2, 10.0, 128).unwrap();
    assert!(
        render_pressure_frames_pcm16(
            &mut renderer,
            &mut stream,
            &CancelGate::new_clock_free(),
            &mut [0.0; 32],
            16,
        )
        .is_err()
    );
    assert_eq!(renderer.frames_rendered(), 0);
    let (_, summary) = stream.finish().unwrap();
    assert_eq!(summary.frames, 0);
}
