use super::*;
use crate::helmholtz::Medium;

const MEAN: [[C64; 2]; 2] = [
    [C64::new(0.4, -0.2), C64::new(-0.1, 0.3)],
    [C64::new(0.2, 0.1), C64::new(0.3, -0.2)],
];
const REAL: [[C64; 2]; 2] = [
    [C64::new(0.7, 0.2), C64::new(-0.3, 0.1)],
    [C64::new(-0.2, 0.4), C64::new(0.1, 0.2)],
];
const PROPER: [[C64; 2]; 2] = [
    [C64::new(-0.1, 0.5), C64::new(0.4, 0.2)],
    [C64::new(0.3, -0.1), C64::new(-0.2, 0.3)],
];

fn exposure_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 4096,
        maximum_receiver_frequency_pairs: 4096,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn realization_limits() -> HarmonicRealizationLimits {
    HarmonicRealizationLimits {
        maximum_lines: 32,
        maximum_field_products: 4096,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn sample_limits() -> HarmonicSampleLimits {
    HarmonicSampleLimits {
        maximum_frames: 8192,
        maximum_receiver_line_samples: 65_536,
    }
}

fn line<'a>(
    owner: &'a PreparedHarmonicObservation,
    coefficients: &'a [C64],
) -> HarmonicExposureLine<'a> {
    HarmonicExposureLine {
        observation: owner,
        coefficients,
    }
}

fn observation(
    omega: f64,
    receivers: &[ReceiverSpec],
    fields: &[[C64; 4]],
) -> PreparedHarmonicObservation {
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let sources: Vec<_> = (0..fields.len() / receivers.len())
        .map(|i| SourceId::PointVolumeFlow(format!("flow-{i}")))
        .collect();
    let mut samples = Vec::new();
    for field in fields {
        samples.push(field[0]);
        samples.extend(
            field[1..]
                .iter()
                .map(|&u| C64::new(0.0, omega * medium.density) * u),
        );
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &sources,
        receivers,
        ContentHash([19; 32]),
        samples,
    )
    .unwrap()
}

// Independent outgoing monopoles: p = -i omega rho Q exp(ikr)/(4 pi r),
// u = p [1/(rho c) + i/(omega rho r)] r_hat. Near-field velocity is retained.
fn monopoles(omega: f64) -> (PreparedHarmonicObservation, [[C64; 4]; 4]) {
    let sources = [[0.0, 0.0, 0.0], [0.25, -0.1, 0.15]];
    let receivers = [
        ReceiverSpec {
            id: "near".into(),
            position_m: [1.0, 0.3, 0.2],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [-0.6, 1.2, 0.7],
        },
    ];
    let mut fields = [[C64::ZERO; 4]; 4];
    for (source, position) in sources.iter().enumerate() {
        for (receiver, spec) in receivers.iter().enumerate() {
            let d: [f64; 3] = core::array::from_fn(|c| spec.position_m[c] - position[c]);
            let r = d.iter().map(|v| v * v).sum::<f64>().sqrt();
            let angle = omega * r / 343.0;
            let green =
                C64::new(angle.cos(), angle.sin()).scale(1.0 / (4.0 * core::f64::consts::PI * r));
            let p = C64::new(0.0, -omega * 1.2) * green;
            let u = p * C64::new(1.0 / (1.2 * 343.0), 1.0 / (omega * 1.2 * r));
            fields[source * 2 + receiver] =
                [p, u.scale(d[0] / r), u.scale(d[1] / r), u.scale(d[2] / r)];
        }
    }
    (observation(omega, &receivers, &fields), fields)
}

fn scene<'a>(
    a: &'a PreparedHarmonicObservation,
    b: &PreparedHarmonicObservation,
) -> PreparedHarmonicExposure<'a> {
    PreparedHarmonicExposure::new_with_statistics(
        a,
        "original-source-clock",
        &[line(a, &MEAN[0]), line(b, &MEAN[1])],
        &[
            &[line(a, &REAL[0]), line(b, &REAL[1])],
            &[line(a, &PROPER[0]), line(b, &PROPER[1])],
        ],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        exposure_limits(),
    )
    .unwrap()
}

fn components(sample: HarmonicFieldSample) -> [f64; 4] {
    [
        sample.pressure_pa,
        sample.particle_velocity_m_s[0],
        sample.particle_velocity_m_s[1],
        sample.particle_velocity_m_s[2],
    ]
}

fn products(v: [f64; 4]) -> [f64; 14] {
    [
        v[0],
        v[1],
        v[2],
        v[3],
        v[0] * v[0],
        v[0] * v[1],
        v[0] * v[2],
        v[0] * v[3],
        v[1] * v[1],
        v[1] * v[2],
        v[1] * v[3],
        v[2] * v[2],
        v[2] * v[3],
        v[3] * v[3],
    ]
}

fn moments(v: HarmonicExposureMoments) -> [f64; 14] {
    [
        v.mean_pressure_pa,
        v.mean_particle_velocity_m_s[0],
        v.mean_particle_velocity_m_s[1],
        v.mean_particle_velocity_m_s[2],
        v.mean_square_pressure_pa2,
        v.mean_intensity_w_m2[0],
        v.mean_intensity_w_m2[1],
        v.mean_intensity_w_m2[2],
        v.velocity_second_m2_s2[0][0],
        v.velocity_second_m2_s2[0][1],
        v.velocity_second_m2_s2[0][2],
        v.velocity_second_m2_s2[1][1],
        v.velocity_second_m2_s2[1][2],
        v.velocity_second_m2_s2[2][2],
    ]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn one_draw_preserves_spherical_wave_phase_sources_and_all_receiver_components() {
    let (a, fa) = monopoles(47.0);
    let (b, fb) = monopoles(73.0);
    let exposure = scene(&a, &b);
    let real = -0.75;
    let proper = C64::new(0.4, -1.2);
    let draws = [
        HarmonicFactorDraw::Real(real),
        HarmonicFactorDraw::ProperComplex(proper),
    ];
    let realized = exposure.realize(&draws, realization_limits()).unwrap();
    drop(exposure);
    assert_eq!(realized.line_count(), 2);
    assert_eq!(realized.factor_draws(), &draws);
    assert_eq!(realized.receivers(), a.receivers());
    assert_eq!(realized.sources(), a.sources());
    assert_eq!(realized.producer_hash(), a.producer_hash());
    assert_eq!(realized.source_clock(), "original-source-clock");
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.017), 8192).unwrap();
    let mut output = vec![HarmonicFieldSample::default(); 66];
    realized
        .sample_fields_into(clock, 0, &mut output, sample_limits())
        .unwrap();
    for (frame, samples) in output.chunks_exact(2).enumerate() {
        let time = 0.017 + frame as f64 / 8192.0;
        for (receiver, sample) in samples.iter().enumerate() {
            let mut expected = [0.0; 4];
            for (frequency, (omega, fields)) in [(47.0, fa), (73.0, fb)].iter().enumerate() {
                let phase = C64::new((omega * time).cos(), -(omega * time).sin());
                for source in 0..2 {
                    let coefficient = MEAN[frequency][source]
                        + REAL[frequency][source].scale(real)
                        + PROPER[frequency][source] * proper;
                    for c in 0..4 {
                        expected[c] += (fields[source * 2 + receiver][c] * coefficient * phase).re;
                    }
                }
            }
            for (actual, expected) in components(*sample).into_iter().zip(expected) {
                close(actual, expected, 3e-13 * (1.0 + expected.abs()));
            }
        }
    }
}

#[test]
fn exact_mixed_factor_ensemble_and_time_integration_recover_exposure_moments() {
    let (a, _) = monopoles(47.0);
    let (b, _) = monopoles(73.0);
    let exposure = scene(&a, &b);
    let start = 0.017;
    let n = 1024usize;
    let clock = HarmonicSampleClock::new(Dd::from_f64(start), 8192).unwrap();
    let window = HarmonicWindow::new(start, n as f64 / 8192.0).unwrap();
    let expected = exposure.observe_window(window).unwrap();
    let mut integral = [[0.0; 14]; 2];
    let mut output = vec![HarmonicFieldSample::default(); 2 * (n + 1)];
    // The eight equally weighted draws exactly realize the mixed second moments.
    for real in [-1.0, 1.0] {
        for proper in [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)] {
            let realized = exposure
                .realize(
                    &[
                        HarmonicFactorDraw::Real(real),
                        HarmonicFactorDraw::ProperComplex(proper),
                    ],
                    realization_limits(),
                )
                .unwrap();
            realized
                .sample_fields_into(clock, 0, &mut output, sample_limits())
                .unwrap();
            for (index, frame) in output.chunks_exact(2).enumerate() {
                let weight = if index == 0 || index == n {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for (sum, sample) in integral.iter_mut().zip(frame) {
                    for (value, moment) in sum.iter_mut().zip(products(components(*sample))) {
                        *value += weight * moment / (24.0 * n as f64);
                    }
                }
            }
        }
    }
    for (actual, expected) in integral.into_iter().zip(expected.points()) {
        for (actual, expected) in actual.into_iter().zip(moments(*expected)) {
            close(actual, expected, 3e-10 * (1.0 + expected.abs()));
        }
    }
}

#[test]
fn absolute_integer_clock_preserves_low_epoch_and_block_partition() {
    let (a, _) = monopoles(47.0);
    let (b, _) = monopoles(73.0);
    let exposure = scene(&a, &b);
    let realized = exposure
        .realize(
            &[
                HarmonicFactorDraw::Real(0.5),
                HarmonicFactorDraw::ProperComplex(C64::new(0.3, -0.8)),
            ],
            realization_limits(),
        )
        .unwrap();
    let first = (1u64 << 54) + 7;
    let clock =
        HarmonicSampleClock::new(Dd::from_pair(-17_592_186_044_416.0, 1.0 / 4096.0), 1024).unwrap();
    assert_eq!(clock.time_at(first).unwrap().to_f64(), 29.0 / 4096.0);
    let mut full = vec![HarmonicFieldSample::default(); 74];
    let mut blocks = full.clone();
    realized
        .sample_fields_into(clock, first, &mut full, sample_limits())
        .unwrap();
    realized
        .sample_fields_into(clock, first, &mut blocks[..10], sample_limits())
        .unwrap();
    realized
        .sample_fields_into(clock, first + 5, &mut blocks[10..38], sample_limits())
        .unwrap();
    realized
        .sample_fields_into(clock, first + 19, &mut blocks[38..], sample_limits())
        .unwrap();
    for (a, b) in full.iter().zip(&blocks) {
        assert_eq!(
            components(*a).map(f64::to_bits),
            components(*b).map(f64::to_bits)
        );
    }
    // This equivalent small epoch is an independent reference for all low bits.
    let local = HarmonicSampleClock::new(Dd::from_f64(29.0 / 4096.0), 1024).unwrap();
    realized
        .sample_fields_into(local, 0, &mut blocks, sample_limits())
        .unwrap();
    for (a, b) in full.iter().zip(&blocks) {
        assert_eq!(
            components(*a).map(f64::to_bits),
            components(*b).map(f64::to_bits)
        );
    }
    assert_ne!(full[0], full[2]);
}

#[test]
fn standing_wave_retains_actual_velocity_and_negative_subperiod_energy_flux() {
    let omega = core::f64::consts::TAU * 16.0;
    let impedance = 1.2 * 343.0;
    let amplitude = 2.0f64.sqrt();
    let receiver = [ReceiverSpec {
        id: "standing-wave".into(),
        position_m: [343.0 / 128.0, 0.0, 0.0],
    }];
    let owner = observation(
        omega,
        &receiver,
        &[[
            C64::from_re(amplitude),
            C64::new(0.0, amplitude / impedance),
            C64::ZERO,
            C64::ZERO,
        ]],
    );
    let exposure = PreparedHarmonicExposure::new(
        &owner,
        "standing",
        &[line(&owner, &[C64::ONE])],
        &[],
        exposure_limits(),
    )
    .unwrap();
    let realized = exposure.realize(&[], realization_limits()).unwrap();
    let clock = HarmonicSampleClock::new(Dd::from_f64(1.0 / 64.0), 65_536).unwrap();
    let n = 1024;
    let mut output = vec![HarmonicFieldSample::default(); n + 1];
    realized
        .sample_fields_into(clock, 0, &mut output, sample_limits())
        .unwrap();
    let mut intensity = 0.0;
    for (i, sample) in output.iter().enumerate() {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        intensity +=
            weight * sample.pressure_pa * sample.particle_velocity_m_s[0] / (3.0 * n as f64);
    }
    close(intensity, -2.0 / (core::f64::consts::PI * impedance), 1e-14);
    assert!(output[n / 2].pressure_pa * output[n / 2].particle_velocity_m_s[0] < 0.0);
    assert!(output[0].pressure_pa.abs() < 1e-14);
    close(
        output[0].particle_velocity_m_s[0],
        amplitude / impedance,
        1e-17,
    );
}

#[test]
fn weak_coherent_residuals_and_distinct_sub_ulp_phase_beats_survive() {
    let receiver = [ReceiverSpec {
        id: "probe".into(),
        position_m: [1.0, 0.0, 0.0],
    }];
    let fields = [[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]];
    let a = observation(1.0, &receiver, &fields);
    let b = observation(1.0 + f64::EPSILON, &receiver, &fields);
    let large = [C64::from_re(1e16)];
    let weak = [C64::ONE];
    let negative = [C64::from_re(-1e16)];
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &a,
        "coherent",
        &[line(&a, &large)],
        &[&[line(&a, &weak)], &[line(&a, &negative)]],
        &[HarmonicExposureFactorStatistics::Real; 2],
        exposure_limits(),
    )
    .unwrap();
    let realized = exposure
        .realize(&[HarmonicFactorDraw::Real(1.0); 2], realization_limits())
        .unwrap();
    assert_eq!(realized.line_count(), 1);
    let mut sample = [HarmonicFieldSample::default()];
    realized
        .sample_fields_into(
            HarmonicSampleClock::new(Dd::ZERO, 8).unwrap(),
            0,
            &mut sample,
            sample_limits(),
        )
        .unwrap();
    assert_eq!(sample[0].pressure_pa, 1.0);

    let beat = PreparedHarmonicExposure::new(
        &a,
        "close-frequencies",
        &[line(&a, &large), line(&b, &negative)],
        &[],
        exposure_limits(),
    )
    .unwrap();
    let realized = beat.realize(&[], realization_limits()).unwrap();
    assert_eq!(realized.line_count(), 2);
    realized
        .sample_fields_into(
            HarmonicSampleClock::new(Dd::ZERO, 8).unwrap(),
            1,
            &mut sample,
            sample_limits(),
        )
        .unwrap();
    let half_beat = 0.5 * f64::EPSILON * 0.125;
    let expected = 2e16 * (0.125 + half_beat).sin() * half_beat.sin();
    close(sample[0].pressure_pa, expected, 2e-15);
}

#[test]
fn sampling_refuses_aliasing_even_at_nulls_and_preserves_output_on_admission_failures() {
    let (owner, _) = monopoles(4000.0);
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owner,
        "bandwidth",
        &[],
        &[&[line(&owner, &MEAN[0])]],
        &[HarmonicExposureFactorStatistics::Real],
        exposure_limits(),
    )
    .unwrap();
    assert!(exposure.realize(&[], realization_limits()).is_err());
    assert!(
        exposure
            .realize(
                &[HarmonicFactorDraw::ProperComplex(C64::ONE)],
                realization_limits()
            )
            .is_err()
    );
    assert!(
        exposure
            .realize(&[HarmonicFactorDraw::Real(f64::NAN)], realization_limits())
            .is_err()
    );
    for limits in [
        HarmonicRealizationLimits {
            maximum_lines: 0,
            ..realization_limits()
        },
        HarmonicRealizationLimits {
            maximum_field_products: 0,
            ..realization_limits()
        },
        HarmonicRealizationLimits {
            maximum_workspace_bytes: 0,
            ..realization_limits()
        },
    ] {
        assert!(
            exposure
                .realize(&[HarmonicFactorDraw::Real(0.0)], limits)
                .is_err()
        );
    }
    assert!(matches!(
        exposure.realize_cancellable(
            &[HarmonicFactorDraw::Real(0.0)],
            realization_limits(),
            || true
        ),
        Err(TransferError::Cancelled)
    ));
    let realized = exposure
        .realize(&[HarmonicFactorDraw::Real(0.0)], realization_limits())
        .unwrap();
    assert_eq!(realized.line_count(), 1);
    let sentinel = HarmonicFieldSample {
        pressure_pa: 123.0,
        particle_velocity_m_s: [4.0, 5.0, 6.0],
    };
    let mut output = [sentinel; 4];
    let aliased = HarmonicSampleClock::new(Dd::ZERO, 1000).unwrap();
    assert!(
        realized
            .sample_fields_into(aliased, 0, &mut output, sample_limits())
            .is_err()
    );
    assert!(realized.validate_sample_range(aliased, 0, 0).is_err());
    assert_eq!(output, [sentinel; 4]);
    let clock = HarmonicSampleClock::new(Dd::ZERO, 4096).unwrap();
    assert!(
        realized
            .sample_fields_into(clock, 0, &mut output[..3], sample_limits())
            .is_err()
    );
    assert!(
        realized
            .sample_fields_into(
                clock,
                0,
                &mut output,
                HarmonicSampleLimits {
                    maximum_frames: 0,
                    ..sample_limits()
                }
            )
            .is_err()
    );
    assert!(
        realized
            .sample_fields_into(
                clock,
                0,
                &mut output,
                HarmonicSampleLimits {
                    maximum_receiver_line_samples: 0,
                    ..sample_limits()
                }
            )
            .is_err()
    );
    assert!(
        realized
            .sample_fields_into(clock, u64::MAX, &mut output, sample_limits())
            .is_err()
    );
    assert!(matches!(
        realized.sample_fields_into_cancellable(clock, 0, &mut output, sample_limits(), || true),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(output, [sentinel; 4]);
    realized
        .sample_fields_into(clock, 0, &mut output, sample_limits())
        .unwrap();
    assert_eq!(output, [HarmonicFieldSample::default(); 4]);
    assert!(HarmonicSampleClock::new(Dd::ZERO, 0).is_err());
    assert!(HarmonicSampleClock::new(Dd::from_f64(f64::INFINITY), 1).is_err());
    assert!(HarmonicSampleClock::new(Dd::from_pair(1e300, 1e280), 48_000).is_err());
}

#[test]
fn empty_exposure_still_retains_physical_scope_and_samples_silence() {
    let (owner, _) = monopoles(47.0);
    let exposure =
        PreparedHarmonicExposure::new(&owner, "silent-clock", &[], &[&[]], exposure_limits())
            .unwrap();
    let realized = exposure
        .realize(
            &[HarmonicFactorDraw::ProperComplex(C64::new(2.0, -3.0))],
            HarmonicRealizationLimits {
                maximum_lines: 0,
                maximum_field_products: 0,
                ..realization_limits()
            },
        )
        .unwrap();
    assert_eq!(realized.line_count(), 0);
    assert_eq!(realized.receivers(), owner.receivers());
    assert_eq!(realized.source_clock(), "silent-clock");
    let mut output = [HarmonicFieldSample {
        pressure_pa: 1.0,
        ..HarmonicFieldSample::default()
    }; 6];
    realized
        .sample_fields_into(
            HarmonicSampleClock::new(Dd::ZERO, 2).unwrap(),
            0,
            &mut output,
            HarmonicSampleLimits {
                maximum_frames: 3,
                maximum_receiver_line_samples: 0,
            },
        )
        .unwrap();
    assert_eq!(output, [HarmonicFieldSample::default(); 6]);
}
