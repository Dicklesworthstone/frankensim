//! Physical waveform members of a shared mean/real/proper source ensemble.
use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    HarmonicEnsembleBudget, HarmonicFactorStatistics, HarmonicFieldFactor, HarmonicFieldLine,
    HarmonicFieldSample, HarmonicSampleBudget, HarmonicSampleClock, PreparedHarmonicEnsemble,
    array::{
        HarmonicArrayBudget, HarmonicArrayReceiver, HarmonicFactorCoordinate, HarmonicFactorDraw,
        HarmonicReadoutBudget, HarmonicReadoutTap, HarmonicRealizationBudget,
        PreparedHarmonicArray,
    },
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::TAU;

const COORDINATES: [HarmonicFactorCoordinate<'static>; 2] = [
    HarmonicFactorCoordinate {
        id: "shared amplitude",
        statistics: HarmonicFactorStatistics::Real,
    },
    HarmonicFactorCoordinate {
        id: "shared proper coefficient",
        statistics: HarmonicFactorStatistics::ProperComplex,
    },
];
const ARRAY: HarmonicArrayBudget = HarmonicArrayBudget {
    max_receivers: 2,
    max_factors: 2,
    max_total_lines: 24,
    max_line_pairs: 256,
};
const REALIZATION: HarmonicRealizationBudget = HarmonicRealizationBudget {
    max_channels: 2,
    max_factors: 2,
    max_total_lines: 24,
    max_line_pairs: 256,
};
const READOUT: HarmonicReadoutBudget = HarmonicReadoutBudget {
    max_taps: 4,
    max_total_lines: 24,
    max_line_pairs: 256,
};
const ZERO: HarmonicFieldSample<3> = HarmonicFieldSample {
    pressure_pa: 0.0,
    particle_velocity_m_s: [0.0; 3],
};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 337.0).unwrap()
}
fn line(
    frequency_hz: f64,
    pressure_pa: C64,
    particle_velocity_m_s: [C64; 3],
) -> HarmonicFieldLine<3> {
    HarmonicFieldLine {
        frequency_hz,
        pressure_pa,
        particle_velocity_m_s,
    }
}
fn ensemble(
    mean: &[HarmonicFieldLine<3>],
    factors: &[HarmonicFieldFactor<'_, 3>],
) -> PreparedHarmonicEnsemble<3> {
    PreparedHarmonicEnsemble::new(
        mean,
        factors,
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: 12,
            max_line_pairs: 144,
        },
    )
    .unwrap()
}
fn factor(
    lines: &[HarmonicFieldLine<3>],
    statistics: HarmonicFactorStatistics,
) -> HarmonicFieldFactor<'_, 3> {
    HarmonicFieldFactor { lines, statistics }
}
fn receiver<'a>(x: f64, ensemble: &'a PreparedHarmonicEnsemble<3>) -> HarmonicArrayReceiver<'a, 3> {
    HarmonicArrayReceiver {
        position_m: [x, 0.0, 0.0],
        acoustic: medium(),
        ensemble,
    }
}
fn sampling(samples: usize) -> HarmonicSampleBudget {
    HarmonicSampleBudget {
        max_samples: samples,
        max_tone_samples: samples * 12,
    }
}
fn components(sample: HarmonicFieldSample<3>) -> [f64; 4] {
    [
        sample.pressure_pa,
        sample.particle_velocity_m_s[0],
        sample.particle_velocity_m_s[1],
        sample.particle_velocity_m_s[2],
    ]
}
fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 2e-10 * scale.max(1e-18),
        "actual {a:e}, reference {b:e}, scale {scale:e}"
    );
}

#[test]
fn shared_realizations_reproduce_spatial_moments_from_physical_waveform_samples() {
    let mean = [
        [line(
            173.0,
            C64::new(0.7, -0.4),
            [
                C64::new(0.002, 0.001),
                C64::new(0.001, -0.003),
                C64::new(-0.003, 0.002),
            ],
        )],
        [line(
            177.0,
            C64::new(-0.2, 0.6),
            [
                C64::new(-0.001, 0.002),
                C64::new(0.002, 0.004),
                C64::new(0.001, -0.002),
            ],
        )],
    ];
    let real = [
        [line(
            181.0,
            C64::new(0.3, 0.2),
            [C64::new(0.004, -0.001); 3],
        )],
        [line(
            173.0,
            C64::new(-0.1, 0.4),
            [C64::new(-0.002, 0.003); 3],
        )],
    ];
    let proper = [
        [line(
            177.0,
            C64::new(0.2, -0.3),
            [C64::new(0.002, 0.004); 3],
        )],
        [line(
            181.0,
            C64::new(0.4, 0.1),
            [C64::new(0.003, -0.001); 3],
        )],
    ];
    let fields: [_; 2] = core::array::from_fn(|r| {
        ensemble(
            &mean[r],
            &[
                factor(&real[r], HarmonicFactorStatistics::Real),
                factor(&proper[r], HarmonicFactorStatistics::ProperComplex),
            ],
        )
    });
    let array = PreparedHarmonicArray::new(
        "physical clock",
        &COORDINATES,
        &[receiver(0.0, &fields[0]), receiver(0.25, &fields[1])],
        ARRAY,
    )
    .unwrap();
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.00013), 250_000).unwrap();
    let window = HarmonicWindow::new(0.00013, 0.004).unwrap();
    let exact = array.observe_pair([0, 1], window).unwrap();
    let mut products = [[0.0; 4]; 4];
    let mut averages = [[0.0; 4]; 2];
    let mut average_products = [[0.0; 4]; 4];
    let mut output = vec![ZERO; 2 * 1001];
    for real_draw in [-1.0, 1.0] {
        for proper_draw in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
            let realized = array
                .realize(
                    &[
                        HarmonicFactorDraw::Real(real_draw),
                        HarmonicFactorDraw::ProperComplex(proper_draw),
                    ],
                    REALIZATION,
                )
                .unwrap();
            realized
                .sample_fields_into(clock, 0, &mut output, sampling(2002))
                .unwrap();
            let mut means = [[0.0; 4]; 2];
            for (index, frame) in output.chunks_exact(2).enumerate() {
                let weight = if index == 0 || index == 1000 {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                } / 3000.0;
                let a = components(frame[0]);
                let b = components(frame[1]);
                for i in 0..4 {
                    means[0][i] += weight * a[i];
                    means[1][i] += weight * b[i];
                    for j in 0..4 {
                        products[i][j] += weight * a[i] * b[j] / 8.0;
                    }
                }
            }
            for i in 0..4 {
                averages[0][i] += means[0][i] / 8.0;
                averages[1][i] += means[1][i] / 8.0;
                for j in 0..4 {
                    average_products[i][j] += means[0][i] * means[1][j] / 8.0;
                }
            }
        }
    }
    let to_matrix = |m: fs_aeroac::observation::array::AcousticCrossMoments<3>| {
        let mut v = [[0.0; 4]; 4];
        v[0][0] = m.pressure_pressure_pa2;
        for i in 0..3 {
            v[0][i + 1] = m.pressure_velocity_w_m2[i];
            v[i + 1][0] = m.velocity_pressure_w_m2[i];
            for j in 0..3 {
                v[i + 1][j + 1] = m.velocity_velocity_m2_s2[i][j];
            }
        }
        v
    };
    let second = to_matrix(exact.second_moments);
    let covariance = to_matrix(exact.averaged_covariance);
    for i in 0..4 {
        for j in 0..4 {
            let scale = if i == 0 && j == 0 {
                1.0
            } else if i == 0 || j == 0 {
                0.01
            } else {
                1e-4
            };
            close(products[i][j], second[i][j], scale);
            close(
                average_products[i][j] - averages[0][i] * averages[1][j],
                covariance[i][j],
                scale,
            );
        }
    }
}

#[test]
fn physical_microphones_reuse_the_same_nonunit_draws_velocity_and_causal_delays() {
    let mean = [line(
        193.0,
        C64::ZERO,
        [C64::new(0.0, 0.004), C64::ZERO, C64::ZERO],
    )];
    let real = [line(
        193.0,
        C64::new(0.2, -0.1),
        [C64::new(-0.001, 0.002); 3],
    )];
    let proper = [line(211.0, C64::new(0.4, 0.3), [C64::new(0.003, 0.001); 3])];
    let a = ensemble(
        &mean,
        &[
            factor(&real, HarmonicFactorStatistics::Real),
            factor(&proper, HarmonicFactorStatistics::ProperComplex),
        ],
    );
    let b = ensemble(
        &proper,
        &[
            factor(&mean, HarmonicFactorStatistics::Real),
            factor(&real, HarmonicFactorStatistics::ProperComplex),
        ],
    );
    let array = PreparedHarmonicArray::new(
        "shared clock",
        &COORDINATES,
        &[receiver(0.0, &a), receiver(0.4, &b)],
        ARRAY,
    )
    .unwrap();
    let taps = [
        HarmonicReadoutTap {
            receiver: 0,
            gain: 1.4,
            delay_s: 0.00027,
            pressure_fraction: 0.0,
            front_axis: [2.0, 0.0, 0.0],
        },
        HarmonicReadoutTap {
            receiver: 1,
            gain: -0.6,
            delay_s: 0.00051,
            pressure_fraction: 0.5,
            front_axis: [0.0, -3.0, 0.0],
        },
    ];
    let readout = array.prepare_readout(&taps, READOUT).unwrap();
    let draws = [
        HarmonicFactorDraw::Real(-2.25),
        HarmonicFactorDraw::ProperComplex(C64::new(0.4, -1.7)),
    ];
    let realized = array.realize(&draws, REALIZATION).unwrap();
    assert_eq!(realized.factor_draws(), draws);
    assert_eq!(realized.receivers()[0].acoustic, medium());
    let channel = realized
        .realize_readout(&readout, REALIZATION, || false)
        .unwrap();
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.001), 8000).unwrap();
    let mut samples = [0.0; 16];
    channel
        .sample_pressure_into(clock, 0, &mut samples, sampling(16))
        .unwrap();
    let z = medium().characteristic_impedance_pa_s_m();
    for (i, &actual) in samples.iter().enumerate() {
        let time = clock.time_at(i as u64).unwrap();
        let a = realized.receivers()[0]
            .field
            .field_at_clock_dd(time - Dd::from_f64(taps[0].delay_s))
            .unwrap();
        let b = realized.receivers()[1]
            .field
            .field_at_clock_dd(time - Dd::from_f64(taps[1].delay_s))
            .unwrap();
        let reference = 1.4 * (-z * a.particle_velocity_m_s[0])
            - 0.6 * (0.5 * b.pressure_pa + 0.5 * z * b.particle_velocity_m_s[1]);
        close(actual, reference, 10.0);
    }
    assert!(samples.iter().any(|x| x.abs() > 0.1));
}

#[test]
fn realization_retains_weak_pressure_velocity_residuals_and_all_original_frequencies() {
    let mean = [
        line(100.0, C64::from_re(1e16), [C64::from_re(1e12); 3]),
        line(100.0, C64::ONE, [C64::from_re(0.25); 3]),
    ];
    let noisy = [line(100.0, C64::from_re(-1e16), [C64::from_re(-1e12); 3])];
    let e = ensemble(&mean, &[factor(&noisy, HarmonicFactorStatistics::Real)]);
    let array = PreparedHarmonicArray::new(
        "cancellation clock",
        &COORDINATES[..1],
        &[receiver(0.0, &e)],
        ARRAY,
    )
    .unwrap();
    let realized = array
        .realize(&[HarmonicFactorDraw::Real(1.0)], REALIZATION)
        .unwrap();
    assert_eq!(
        realized.receivers()[0]
            .field
            .field_at_clock_dd(Dd::ZERO)
            .unwrap(),
        HarmonicFieldSample {
            pressure_pa: 1.0,
            particle_velocity_m_s: [0.25; 3],
        }
    );
    let mean = [line(100.0, C64::ONE, [C64::ZERO; 3])];
    let high = [line(5000.0, C64::ONE, [C64::ZERO; 3])];
    let anti = [line(5000.0, -C64::ONE, [C64::ZERO; 3])];
    let coordinates = [
        COORDINATES[0],
        HarmonicFactorCoordinate {
            id: "second real coefficient",
            statistics: HarmonicFactorStatistics::Real,
        },
    ];
    let e = ensemble(
        &mean,
        &[
            factor(&high, HarmonicFactorStatistics::Real),
            factor(&anti, HarmonicFactorStatistics::Real),
        ],
    );
    let array =
        PreparedHarmonicArray::new("bandwidth clock", &coordinates, &[receiver(0.0, &e)], ARRAY)
            .unwrap();
    let silent_readout = array
        .prepare_readout(&[HarmonicReadoutTap::pressure(0, 0.0, 0.0)], READOUT)
        .unwrap();
    for draw in [0.0, 1.0] {
        let realized = array
            .realize(&[HarmonicFactorDraw::Real(draw); 2], REALIZATION)
            .unwrap();
        assert_eq!(realized.receivers()[0].field.frequency_count(), 2);
        let sentinel = HarmonicFieldSample {
            pressure_pa: 123.0,
            particle_velocity_m_s: [456.0; 3],
        };
        let mut out = [sentinel; 4];
        let clock = HarmonicSampleClock::new(Dd::ZERO, 8000).unwrap();
        assert!(
            realized
                .sample_fields_into(clock, 0, &mut out, sampling(4))
                .is_err()
        );
        assert_eq!(out, [sentinel; 4]);
        let silent = realized
            .realize_readout(&silent_readout, REALIZATION, || false)
            .unwrap();
        let mut output = [123.0; 4];
        assert!(
            silent
                .sample_pressure_into(clock, 0, &mut output, sampling(4))
                .is_err()
        );
        assert_eq!(output, [123.0; 4]);
    }
}

#[test]
fn integer_clock_blocks_preserve_complex_phase_and_pressure_velocity_sample_identity() {
    let tone = [line(
        128.0,
        C64::new(0.4, 0.7),
        [C64::new(0.003, -0.002); 3],
    )];
    let e = ensemble(
        &[],
        &[factor(&tone, HarmonicFactorStatistics::ProperComplex)],
    );
    let array = PreparedHarmonicArray::new(
        "long clock",
        &COORDINATES[1..],
        &[receiver(0.0, &e), receiver(0.2, &e)],
        ARRAY,
    )
    .unwrap();
    let draw = C64::new(-0.3, 1.4);
    let realized = array
        .realize(&[HarmonicFactorDraw::ProperComplex(draw)], REALIZATION)
        .unwrap();
    let clock = HarmonicSampleClock::new(Dd::from_pair(1_099_511_627_776.0, 0.0625), 4096).unwrap();
    let first = (1_u64 << 53) + 5;
    let mut whole = [ZERO; 64];
    let mut chunked = [ZERO; 64];
    realized
        .sample_fields_into(clock, first, &mut whole, sampling(64))
        .unwrap();
    let mut start = 0;
    for count in [1, 5, 2, 24] {
        realized
            .sample_fields_into(
                clock,
                first + start as u64,
                &mut chunked[2 * start..2 * (start + count)],
                sampling(2 * count),
            )
            .unwrap();
        start += count;
    }
    assert_eq!(whole, chunked);
    let phase = C64::new((TAU * 5.0 / 32.0).cos(), -(TAU * 5.0 / 32.0).sin());
    close(
        whole[0].pressure_pa,
        (tone[0].pressure_pa * draw * phase).re,
        1.0,
    );
    close(
        whole[0].particle_velocity_m_s[0],
        (tone[0].particle_velocity_m_s[0] * draw * phase).re,
        0.01,
    );
    let mut scalar = [0.0; 32];
    let mut field = [ZERO; 32];
    realized.receivers()[0]
        .field
        .sample_pressure_into(clock, first, &mut scalar, sampling(32))
        .unwrap();
    realized.receivers()[0]
        .field
        .sample_field_into(clock, first, &mut field, sampling(32))
        .unwrap();
    for i in 0..32 {
        assert_eq!(scalar[i], whole[2 * i].pressure_pa);
        assert_eq!(field[i], whole[2 * i]);
        assert_eq!(whole[2 * i], whole[2 * i + 1]);
    }
}

#[test]
fn realization_binding_budgets_clock_and_cancellation_are_explicit() {
    let tone = [line(137.0, C64::ONE, [C64::new(0.002, 0.001); 3])];
    let e = ensemble(&tone, &[factor(&tone, HarmonicFactorStatistics::Real)]);
    let receivers = [receiver(0.0, &e), receiver(0.2, &e)];
    let array = PreparedHarmonicArray::new("clock", &COORDINATES[..1], &receivers, ARRAY).unwrap();
    for draws in [
        &[][..],
        &[HarmonicFactorDraw::ProperComplex(C64::ONE)][..],
        &[HarmonicFactorDraw::Real(f64::NAN)][..],
    ] {
        assert!(array.realize(draws, REALIZATION).is_err());
    }
    let draws = [HarmonicFactorDraw::Real(0.5)];
    assert!(
        array
            .realize(
                &draws,
                HarmonicRealizationBudget {
                    max_line_pairs: 7,
                    ..REALIZATION
                }
            )
            .is_err()
    );
    assert!(matches!(
        array.realize_cancellable(&draws, REALIZATION, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let realized = array.realize(&draws, REALIZATION).unwrap();
    let other_clock =
        PreparedHarmonicArray::new("other clock", &COORDINATES[..1], &receivers, ARRAY).unwrap();
    let other_id = [HarmonicFactorCoordinate {
        id: "unrelated real coefficient",
        statistics: HarmonicFactorStatistics::Real,
    }];
    let other_source = PreparedHarmonicArray::new("clock", &other_id, &receivers, ARRAY).unwrap();
    for a in [&other_clock, &other_source] {
        let channel = a
            .prepare_readout(&[HarmonicReadoutTap::pressure(0, 1.0, 0.0)], READOUT)
            .unwrap();
        assert!(
            realized
                .realize_readout(&channel, REALIZATION, || false)
                .is_err()
        );
    }
    let clock = HarmonicSampleClock::new(Dd::ZERO, 8000).unwrap();
    let mut wrong_shape = [ZERO; 3];
    assert!(
        realized
            .sample_fields_into(clock, 0, &mut wrong_shape, sampling(3))
            .is_err()
    );
    let mut out = [ZERO; 32];
    assert!(
        realized
            .sample_fields_into(clock, u64::MAX, &mut out, sampling(32))
            .is_err()
    );
    assert!(
        realized
            .sample_fields_into(
                clock,
                0,
                &mut out,
                HarmonicSampleBudget {
                    max_samples: 32,
                    max_tone_samples: 31
                }
            )
            .is_err()
    );
    assert_eq!(out, [ZERO; 32]);
    let mut polls = 0;
    assert!(matches!(
        realized.sample_fields_into_cancellable(clock, 0, &mut out, sampling(32), || {
            polls += 1;
            polls == 25
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    realized
        .sample_fields_into(clock, 0, &mut out, sampling(32))
        .unwrap();
    let mut replay = [ZERO; 32];
    realized
        .sample_fields_into(clock, 0, &mut replay, sampling(32))
        .unwrap();
    assert_eq!(out, replay);
}
