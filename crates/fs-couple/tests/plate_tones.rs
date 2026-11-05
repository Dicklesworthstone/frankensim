//! Continuous plate forces reach the existing calibrated pressure/PCM consumer.
use fs_couple::pcm_wav::{
    encode_pcm16_wav,
    observation::PressureRenderer,
    stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16},
};
use fs_couple::render::{
    RenderError,
    plate::tones::{VkToneEnvelope, VkToneGainKnot, VkTonePressureRenderer, VkToneRenderConfig},
};
use fs_couple::thin_plate::{PreparedVkTones, VkBody, VkToneConfig};
use fs_exec::CancelGate;
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_phs::{ResistanceWitness, linear_exact::DenseExponentialBudget};
use std::io::Cursor;

const RATE: u32 = 1_000;
const SAMPLES: usize = 257;
const OMEGAS: [f64; 2] = [11.0, 17.0];
const COSINE: [f64; 2] = [4.0, -1.2];
const SINE: [f64; 2] = [1.3, 2.0];
const DAMPING: f64 = 0.4;
const KNOTS: [VkToneGainKnot; 6] = [
    VkToneGainKnot {
        sample: 0,
        gain: 0.0,
    },
    VkToneGainKnot {
        sample: 47,
        gain: 1.0,
    },
    VkToneGainKnot {
        sample: 95,
        gain: 0.0,
    },
    VkToneGainKnot {
        sample: 173,
        gain: 0.0,
    },
    VkToneGainKnot {
        sample: 221,
        gain: -0.7,
    },
    VkToneGainKnot {
        sample: 257,
        gain: 0.0,
    },
];

fn source(nonlinear: bool, affine: bool) -> PreparedVkTones {
    let body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![8.0],
            channels: if nonlinear {
                vec![StressChannel {
                    coefficient: 300.0,
                    coupling: vec![1.0],
                }]
            } else {
                vec![]
            },
        },
        &[DAMPING],
        &[1.0],
        vec![0.02],
    )
    .unwrap();
    let channel = ResistanceWitness::try_from_source_rows(
        2,
        vec![0.0, 0.0, 0.0, DAMPING],
        vec![0.0, DAMPING.sqrt()],
    )
    .unwrap();
    let config = VkToneConfig {
        time_step_s: f64::from(RATE).recip(),
        nyquist_guard_fraction: 0.9,
        dense_budget: DenseExponentialBudget {
            max_dimension: 128,
            max_dense_bytes: 16 * 1024 * 1024,
            max_multiply_terms: 1_000_000_000,
        },
    };
    if affine {
        body.into_prepared_affine_tones(&OMEGAS, config, &[channel])
            .unwrap()
    } else {
        body.into_prepared_tones(&OMEGAS, config, &[channel])
            .unwrap()
    }
}

fn config() -> VkToneRenderConfig {
    VkToneRenderConfig {
        sample_rate_hz: RATE,
        samples: SAMPLES as u64,
        max_block_len: 128,
        max_driver_bytes: 4096,
        max_envelope_knots: 32,
        density_kg_m3: 1.2,
        distance_m: 1.7,
    }
}

fn envelopes() -> [VkToneEnvelope<'static>; 2] {
    [
        VkToneEnvelope {
            port: 0,
            tone: 0,
            knots: &KNOTS,
        },
        VkToneEnvelope {
            port: 0,
            tone: 1,
            knots: &KNOTS,
        },
    ]
}

fn fixed(nonlinear: bool) -> VkTonePressureRenderer<'static> {
    VkTonePressureRenderer::fixed(source(nonlinear, false), &COSINE, &SINE, config()).unwrap()
}

#[test]
fn fixed_forces_preserve_original_state_work_and_heat_across_callback_partitions() {
    for nonlinear in [false, true] {
        let mut physical = source(nonlinear, false);
        let (mut cosine, mut sine) = (COSINE, SINE);
        let mut expected = Vec::new();
        let mut last = None;
        for _ in 0..SAMPLES {
            let frame = physical
                .drive_and_radiate(&cosine, &sine, 1.2, 1.7)
                .unwrap();
            expected.push(frame.pressure_pa);
            last = Some(frame);
            for tone in 0..OMEGAS.len() {
                (cosine[tone], sine[tone]) = physical
                    .advance_force_phase(tone, cosine[tone], sine[tone])
                    .unwrap();
            }
        }
        for block in [1, 7, 31, 128] {
            let mut renderer = fixed(nonlinear);
            let mut pressure = vec![0.0; SAMPLES];
            for chunk in pressure.chunks_mut(block) {
                renderer.block(chunk).unwrap();
            }
            assert_eq!(pressure, expected);
            assert_eq!(renderer.samples_rendered(), SAMPLES as u64);
            assert_eq!(renderer.source().state(), physical.state());
            assert_eq!(renderer.last_step(), last.as_ref());
            assert_eq!(renderer.channel_heat_j(), physical.channel_heat_j());
            assert!(renderer.channel_heat_j()[0] > 0.0);
        }
    }
}

// Independent continuous ODE reference: it evaluates the authored phase from
// absolute physical time, never the renderer's cached rotation or endpoint API.
fn gain(sample: f64) -> f64 {
    for pair in KNOTS.windows(2) {
        if sample <= pair[1].sample as f64 {
            let t = (sample - pair[0].sample as f64) / (pair[1].sample - pair[0].sample) as f64;
            return pair[0].gain * (1.0 - t) + pair[1].gain * t;
        }
    }
    KNOTS.last().unwrap().gain
}

fn derivative(x: [f64; 4], time: f64, nonlinear: bool) -> [f64; 4] {
    let force = gain(time * f64::from(RATE))
        * (0..2)
            .map(|tone| {
                COSINE[tone] * (OMEGAS[tone] * time).cos()
                    + SINE[tone] * (OMEGAS[tone] * time).sin()
            })
            .sum::<f64>();
    [
        x[1],
        force - 64.0 * x[0] - DAMPING * x[1] - if nonlinear { 300.0 * x[0].powi(3) } else { 0.0 },
        force * x[1],
        DAMPING * x[1] * x[1],
    ]
}

fn ode_oracle(nonlinear: bool) -> (Vec<f64>, [f64; 4]) {
    let mut x = [0.0; 4];
    let mut pressures = Vec::new();
    let h = f64::from(RATE).recip() / 16.0;
    for sample in 0..SAMPLES {
        let before = x[1];
        x[2] = 0.0;
        x[3] = 0.0;
        for sub in 0..16 {
            let time = (sample * 16 + sub) as f64 * h;
            let add = |v: [f64; 4], scale: f64| core::array::from_fn(|i| x[i] + scale * v[i]);
            let a = derivative(x, time, nonlinear);
            let b = derivative(add(a, h / 2.0), time + h / 2.0, nonlinear);
            let c = derivative(add(b, h / 2.0), time + h / 2.0, nonlinear);
            let d = derivative(add(c, h), time + h, nonlinear);
            for i in 0..4 {
                x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
            }
        }
        pressures
            .push(1.2 * 0.02 * (x[1] - before) * f64::from(RATE) / (core::f64::consts::TAU * 1.7));
    }
    (pressures, x)
}

#[test]
fn affine_physical_gains_match_continuous_ode_and_preserve_phase_through_silence() {
    for nonlinear in [false, true] {
        let (oracle_pressure, oracle) = ode_oracle(nonlinear);
        let mut previous = None;
        for block in [1, 13, 128] {
            let mut renderer = VkTonePressureRenderer::affine(
                source(nonlinear, true),
                &COSINE,
                &SINE,
                &envelopes(),
                config(),
            )
            .unwrap();
            let mut pressure = vec![0.0; SAMPLES];
            for chunk in pressure.chunks_mut(block) {
                renderer.block(chunk).unwrap();
            }
            let tolerance = if nonlinear { 2e-7 } else { 2e-11 };
            for (actual, expected) in pressure.iter().zip(&oracle_pressure) {
                assert!(
                    (actual - expected).abs() < tolerance,
                    "nonlinear={nonlinear}: {actual} != {expected}"
                );
            }
            for (actual, expected) in renderer.source().state().iter().zip(&oracle) {
                assert!((actual - expected).abs() < tolerance);
            }
            let frame = renderer.last_step().unwrap();
            assert!((frame.work.supplied - oracle[2]).abs() < tolerance);
            assert!((renderer.channel_heat_j()[0] - oracle[3]).abs() < tolerance);
            // A silent force interval still carries the vibrating physical body.
            assert!(pressure[100..170].iter().any(|p| p.abs() > 1e-5));
            if let Some(previous) = &previous {
                assert_eq!(&pressure, previous);
            }
            previous = Some(pressure);
        }
    }
}

#[test]
fn calibrated_pcm_stream_reuses_source_clock_and_cancelled_prefix_without_normalization() {
    let mut reference = fixed(true);
    let mut pressure = vec![0.0; SAMPLES];
    for chunk in pressure.chunks_mut(128) {
        reference.block(chunk).unwrap();
    }
    let (expected, clips) = encode_pcm16_wav(&pressure, RATE, 0.003).unwrap();
    assert!(clips > 0);
    for block in [7, 128] {
        let mut renderer = fixed(true);
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), RATE, 0.003, 128).unwrap();
        let mut scratch = vec![0.0; block];
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert_eq!(
            render_pressure_pcm16(
                &mut renderer,
                &mut stream,
                &cancelled,
                &mut scratch,
                SAMPLES as u64
            )
            .unwrap(),
            ScheduledWavProgress::Cancelled { samples: 0 }
        );
        assert_eq!(renderer.samples_rendered(), 0);
        let gate = CancelGate::new_clock_free();
        for count in [31, SAMPLES as u64 - 31] {
            assert_eq!(
                render_pressure_pcm16(&mut renderer, &mut stream, &gate, &mut scratch, count)
                    .unwrap(),
                ScheduledWavProgress::Completed { samples: count }
            );
        }
        let (sink, summary) = stream.finish().unwrap();
        assert_eq!(sink.into_inner(), expected);
        assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(summary.full_scale_pa, 0.003);
        assert_eq!(renderer.source().state(), reference.source().state());
    }
}

#[test]
fn request_admission_is_read_only_and_observation_failure_poison_preserves_accepted_prefix() {
    let mut renderer = fixed(false);
    let before = renderer.source().state().to_vec();
    assert!(renderer.validate_sample_rate(RATE + 1).is_err());
    assert!(renderer.validate_sample_count(SAMPLES as u64 + 1).is_err());
    let mut too_large = [9.0; 129];
    assert!(renderer.block(&mut too_large).is_err());
    assert_eq!(too_large, [9.0; 129]);
    assert_eq!(renderer.samples_rendered(), 0);
    assert_eq!(renderer.source().state(), before);
    renderer.block(&mut [0.0; 1]).unwrap();

    let knots = [
        VkToneGainKnot {
            sample: 0,
            gain: 0.0,
        },
        VkToneGainKnot {
            sample: 1,
            gain: 0.0,
        },
        VkToneGainKnot {
            sample: 2,
            gain: 1.0,
        },
    ];
    let envelopes = [
        VkToneEnvelope {
            port: 0,
            tone: 0,
            knots: &knots,
        },
        VkToneEnvelope {
            port: 0,
            tone: 1,
            knots: &knots,
        },
    ];
    let mut observer = config();
    observer.density_kg_m3 = 1e300;
    observer.distance_m = 1e-300;
    let mut renderer =
        VkTonePressureRenderer::affine(source(false, true), &COSINE, &SINE, &envelopes, observer)
            .unwrap();
    let mut output = [9.0; 3];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [0.0, 9.0, 9.0]);
    assert_eq!(renderer.samples_rendered(), 1);
    assert_eq!(renderer.source().state(), &[0.0, 0.0]);
    assert_eq!(renderer.last_step().unwrap().work.supplied, 0.0);
    assert_eq!(renderer.channel_heat_j(), &[0.0]);
    assert!(matches!(
        renderer.block(&mut [0.0; 1]),
        Err(RenderError::Poisoned)
    ));
}

#[test]
fn malformed_or_unbudgeted_force_programs_refuse_before_rendering() {
    let mut limited = config();
    limited.max_driver_bytes = 1;
    assert!(VkTonePressureRenderer::fixed(source(false, false), &COSINE, &SINE, limited).is_err());
    limited = config();
    limited.max_envelope_knots = 1;
    assert!(
        VkTonePressureRenderer::affine(source(false, true), &COSINE, &SINE, &envelopes(), limited)
            .is_err()
    );
    let mut invalid = envelopes();
    invalid[1].tone = 0;
    assert!(
        VkTonePressureRenderer::affine(source(false, true), &COSINE, &SINE, &invalid, config())
            .is_err()
    );
    assert!(
        VkTonePressureRenderer::fixed(source(false, false), &[f64::NAN, 0.0], &SINE, config())
            .is_err()
    );
    assert!(VkTonePressureRenderer::fixed(source(false, true), &COSINE, &SINE, config()).is_err());
}
