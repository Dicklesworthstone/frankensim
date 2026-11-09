use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn source_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn limits() -> HarmonicSurfaceFieldLimits {
    HarmonicSurfaceFieldLimits {
        maximum_quadrature_points: 32,
        maximum_temporal_projection_terms: 100_000,
        maximum_quadrature_frequency_pairs: 8192,
    }
}
fn owner(omega: f64, fields: Vec<[C64; 4]>) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 2.0,
            sound_speed: 4.0,
        },
        sources: vec![SourceId::StructuralMode("original-source".into())],
        receivers: fields
            .iter()
            .enumerate()
            .map(|(i, _)| ReceiverSpec {
                id: format!("physical-{i}"),
                position_m: [0.2 * i as f64, 3.0, -0.4],
            })
            .collect(),
        producer_hash: ContentHash([73; 32]),
        transfers: fields,
        radiation_momentum: true,
    }
}
fn bank<'a>(
    owners: &'a [PreparedHarmonicObservation; 3],
    coefficients: [C64; 3],
    proper: bool,
) -> PreparedHarmonicExposure<'a> {
    let coefficients = coefficients.map(|x| [x]);
    let lines: [_; 3] = core::array::from_fn(|i| HarmonicExposureLine {
        observation: &owners[i],
        coefficients: &coefficients[i],
    });
    if proper {
        PreparedHarmonicExposure::new(
            &owners[0],
            "original-clock",
            &[],
            &[&lines],
            source_limits(),
        )
        .unwrap()
    } else {
        PreparedHarmonicExposure::new(&owners[0], "original-clock", &lines, &[], source_limits())
            .unwrap()
    }
}
fn plane(frequencies: [f64; 3]) -> [PreparedHarmonicObservation; 3] {
    frequencies.map(|omega| {
        owner(
            omega,
            vec![[C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO]],
        )
    })
}
fn surface() -> [HarmonicSurfaceQuadrature; 1] {
    [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [2.0, 0.0, 0.0],
    }]
}
fn close(actual: Dd, expected: f64, relative: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= relative * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:e}"
    );
}
fn unpack(x: HarmonicSurfaceFieldMoments) -> [Dd; 7] {
    [
        x.mean_power_w,
        x.mean_linear_momentum_rate_n[0],
        x.mean_linear_momentum_rate_n[1],
        x.mean_linear_momentum_rate_n[2],
        x.mean_angular_momentum_rate_nm[0],
        x.mean_angular_momentum_rate_nm[1],
        x.mean_angular_momentum_rate_nm[2],
    ]
}

#[test]
fn dark_plane_power_stress_and_angular_flux_retain_original_interference() {
    let carrier = 2.0_f64.powi(49);
    let owners = plane([carrier, carrier + 1.0, carrier + 2.0]);
    let exposure = bank(&owners, [1.0, -2.0, 1.0].map(C64::from_re), true);
    let quadrature = surface();
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let start = -0.5 * duration;
        let window = HarmonicWindow::new(start, duration).unwrap();
        let actual = exposure
            .surface_field_window(&quadrature, [0.0; 3], window, limits(), || false)
            .unwrap();
        let square = duration.powi(4) / 160.0;
        close(actual.mean_power_w, 16.0 * square, 3e-12);
        close(actual.mean_linear_momentum_rate_n[0], 4.0 * square, 3e-12);
        close(
            actual.mean_angular_momentum_rate_nm[1],
            -1.6 * square,
            3e-12,
        );
        close(
            actual.mean_angular_momentum_rate_nm[2],
            -12.0 * square,
            3e-12,
        );
        assert_eq!(actual.mean_linear_momentum_rate_n[1..], [Dd::ZERO; 2]);
        let interval = exposure
            .surface_field_interval_from_clock_dd(
                &quadrature,
                [0.0; 3],
                Dd::from_f64(start),
                duration,
                limits(),
                || false,
            )
            .unwrap();
        close(interval.energy_j, 16.0 * square * duration, 3e-12);
        close(
            interval.wave_momentum.linear_momentum_ns[0],
            4.0 * square * duration,
            3e-12,
        );
    }
}

const COEFFICIENTS: [[[f64; 2]; 3]; 3] = [
    [[0.7, 0.2], [-0.4, 0.3], [0.1, -0.6]],
    [[-0.2, 0.1], [0.3, 0.7], [0.6, -0.2]],
    [[0.4, -0.3], [0.8, 0.1], [-0.2, 0.5]],
];

fn oracle(
    owners: &[PreparedHarmonicObservation; 3],
    surface: &[HarmonicSurfaceQuadrature],
    origin: [f64; 3],
    window: HarmonicWindow,
) -> [f64; 7] {
    let n = 2048;
    let mut result = [0.0; 7];
    for sample_index in 0..=n {
        let fraction = sample_index as f64 / n as f64;
        let time = window.start_s() + window.duration_s() * fraction;
        let weight = if sample_index == 0 || sample_index == n {
            1.0
        } else if sample_index % 2 == 0 {
            2.0
        } else {
            4.0
        } * if window.profile() == HarmonicWindowProfile::Triangular {
            2.0 * (1.0 - (2.0 * fraction - 1.0).abs())
        } else {
            1.0
        } / (3.0 * n as f64);
        for sample in surface {
            for (bank, coefficients) in COEFFICIENTS.iter().enumerate() {
                let mut fields = [C64::ZERO; 4];
                for (i, owner) in owners.iter().enumerate() {
                    let phase = C64::new((owner.omega * time).cos(), -(owner.omega * time).sin());
                    let amplitude = C64::new(coefficients[i][0], coefficients[i][1]) * phase;
                    for c in 0..4 {
                        fields[c] = fields[c] + owner.transfers[sample.receiver][c] * amplitude;
                    }
                }
                let moment = |i: usize, j: usize| {
                    if bank == 2 {
                        0.5 * (fields[i].re * fields[j].re + fields[i].im * fields[j].im)
                    } else {
                        fields[i].re * fields[j].re
                    }
                };
                let rho = owners[0].medium.density;
                let c = owners[0].medium.sound_speed;
                let isotropic = 0.5
                    * (moment(0, 0) / (rho * c * c)
                        - rho * (moment(1, 1) + moment(2, 2) + moment(3, 3)));
                let mut traction = [0.0; 3];
                for i in 0..3 {
                    result[0] += weight * sample.oriented_area_m2[i] * moment(0, i + 1);
                    traction[i] = isotropic * sample.oriented_area_m2[i];
                    for j in 0..3 {
                        traction[i] += rho * sample.oriented_area_m2[j] * moment(i + 1, j + 1);
                    }
                    result[i + 1] += weight * traction[i];
                }
                let arm: [f64; 3] = core::array::from_fn(|i| {
                    owners[0].receivers[sample.receiver].position_m[i] - origin[i]
                });
                for i in 0..3 {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    result[i + 4] += weight * (arm[j] * traction[k] - arm[k] * traction[j]);
                }
            }
        }
    }
    result
}

#[test]
fn mixed_real_proper_wave_stress_matches_independent_time_quadrature_and_fallback() {
    let frequencies = [1.2, 1.7, 2.3];
    let owners: [_; 3] = core::array::from_fn(|i| {
        owner(
            frequencies[i],
            (0..2)
                .map(|r| {
                    core::array::from_fn(|c| {
                        C64::new(
                            0.2 + 0.3 * c as f64 - 0.1 * i as f64,
                            0.15 * (1 + r + i) as f64 - 0.12 * c as f64,
                        )
                    })
                })
                .collect(),
        )
    });
    let coefficients = COEFFICIENTS.map(|bank| bank.map(|[r, i]| [C64::new(r, i)]));
    let banks: [[_; 3]; 3] = core::array::from_fn(|b| {
        core::array::from_fn(|i| HarmonicExposureLine {
            observation: &owners[i],
            coefficients: &coefficients[b][i],
        })
    });
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "mixed-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        source_limits(),
    )
    .unwrap();
    let surface = [
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.7, -0.3, 0.2],
        },
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [-0.2, 0.6, 0.1],
        },
    ];
    let origin = [0.3, 0.2, -0.1];
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.2, 0.3, profile).unwrap();
        let actual = unpack(
            exposure
                .surface_field_window(&surface, origin, window, limits(), || false)
                .unwrap(),
        );
        let expected = oracle(&owners, &surface, origin, window);
        for i in 0..7 {
            close(actual[i], expected[i], 3e-12);
        }
        close(
            exposure
                .surface_field_power_window(&surface, window, limits(), || false)
                .unwrap(),
            expected[0],
            3e-12,
        );
    }
}

#[test]
fn dd_epoch_and_partitioned_signed_transport_retain_source_continuation() {
    let owners = plane([1.0, 2.0, 3.0]);
    let exposure = bank(
        &owners,
        [C64::new(0.7, 0.2), C64::new(-0.3, 0.4), C64::new(0.2, -0.6)],
        false,
    );
    let epoch = Dd::from_pair(2.0_f64.powi(49), 0.001);
    let query = |start, dt| {
        exposure
            .surface_field_interval_from_clock_dd(&surface(), [0.0; 3], start, dt, limits(), || {
                false
            })
            .unwrap()
    };
    let whole = query(epoch, 0.04);
    let first = query(epoch, 0.02);
    let last = query(epoch + Dd::from_f64(0.02), 0.02);
    close(
        first.energy_j + last.energy_j,
        whole.energy_j.to_f64(),
        3e-14,
    );
    for i in 0..3 {
        close(
            first.wave_momentum.linear_momentum_ns[i] + last.wave_momentum.linear_momentum_ns[i],
            whole.wave_momentum.linear_momentum_ns[i].to_f64(),
            3e-14,
        );
    }
    let rounded = query(Dd::from_f64(epoch.to_f64()), 0.04);
    assert!((rounded.energy_j - whole.energy_j).to_f64().abs() > 1e-6);
    assert_eq!(epoch.to_f64(), (epoch + Dd::from_f64(0.02)).to_f64());
    close(
        exposure
            .surface_field_energy_from_clock_dd(&surface(), epoch, 0.04, limits(), || false)
            .unwrap(),
        whole.energy_j.to_f64(),
        3e-14,
    );
}

#[test]
fn oriented_area_origin_and_original_medium_authority_are_preserved() {
    let mut owners = plane([1.0, 2.0, 3.0]);
    let exposure = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], true);
    let window = HarmonicWindow::new(0.2, 0.1).unwrap();
    let quadrature = surface();
    let base = exposure
        .surface_field_window(&quadrature, [0.0; 3], window, limits(), || false)
        .unwrap();
    let split = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }; 2];
    let split_result = exposure
        .surface_field_window(&split, [0.0; 3], window, limits(), || false)
        .unwrap();
    // Subdivision changes only the compensated reduction order, not area or
    // physical flux. It need not preserve the last double-double bit.
    for (a, b) in unpack(base).into_iter().zip(unpack(split_result)) {
        assert!((a - b).to_f64().abs() <= 1e-29 * a.to_f64().abs().max(1.0));
    }
    let reverse = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [-2.0, 0.0, 0.0],
    }];
    let reversed = exposure
        .surface_field_window(&reverse, [0.0; 3], window, limits(), || false)
        .unwrap();
    assert_eq!(unpack(base).map(|x| -x), unpack(reversed));
    let shifted = exposure
        .surface_field_window(&quadrature, [0.0, 1.0, 0.0], window, limits(), || false)
        .unwrap();
    assert_eq!(shifted.mean_power_w, base.mean_power_w);
    close(
        shifted.mean_angular_momentum_rate_nm[2],
        (base.mean_angular_momentum_rate_nm[2] + base.mean_linear_momentum_rate_n[0]).to_f64(),
        1e-14,
    );
    owners[2].radiation_momentum = false;
    let relaxing = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], true);
    assert!(
        relaxing
            .surface_field_window(&quadrature, [0.0; 3], window, limits(), || false)
            .is_err()
    );
    assert_eq!(
        relaxing
            .surface_field_power_window(&quadrature, window, limits(), || false)
            .unwrap(),
        base.mean_power_w
    );
}

#[test]
fn interval_scale_precedes_large_products_and_preserves_finite_physical_half() {
    let mut owners = [1.0, 2.0, 3.0].map(|omega| {
        owner(
            omega,
            vec![[
                C64::from_re(1.5e200),
                C64::from_re(1.5e200),
                C64::ZERO,
                C64::ZERO,
            ]],
        )
    });
    for owner in &mut owners {
        owner.medium.density = 1.0;
        owner.medium.sound_speed = 1.0;
    }
    let exposure = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], true);
    let surface = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    let result = exposure
        .surface_field_interval_from_clock_dd(
            &surface,
            [0.0, 3.0, -0.4],
            Dd::ZERO,
            1e-100,
            limits(),
            || false,
        )
        .unwrap();
    close(result.energy_j, 1.125e300, 4e-15);
    close(result.wave_momentum.linear_momentum_ns[0], 1.125e300, 4e-15);
    assert_eq!(result.wave_momentum.angular_momentum_nms, [Dd::ZERO; 3]);
    // A lost proper half-factor must refuse, not fabricate zero energy for
    // these large fields over the smallest representable positive interval.
    assert!(matches!(
        exposure.surface_field_interval_from_clock_dd(
            &surface,
            [0.0, 3.0, -0.4],
            Dd::ZERO,
            f64::from_bits(1),
            limits(),
            || false,
        ),
        Err(TransferError::Invalid(
            "acoustic surface physical time scale is unrepresentable"
        ))
    ));
    assert!(
        exposure
            .surface_field_power_window(
                &surface,
                HarmonicWindow::new(0.0, 1e-100).unwrap(),
                limits(),
                || false
            )
            .is_err()
    );
}

#[test]
fn work_caps_invalid_geometry_and_cancellation_refuse_without_changing_owner() {
    let owners = plane([1.0, 2.0, 3.0]);
    let exposure = bank(&owners, [C64::ONE; 3], false);
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let projection_terms = mode_count(&exposure.mean, window).unwrap().unwrap() * (3 + 4 * 3 + 10);
    let exact = HarmonicSurfaceFieldLimits {
        maximum_quadrature_points: 1,
        maximum_temporal_projection_terms: projection_terms,
        maximum_quadrature_frequency_pairs: 0,
    };
    let query = |caps| exposure.surface_field_window(&surface(), [0.0; 3], window, caps, || false);
    let accepted = query(exact).unwrap();
    assert!(
        query(HarmonicSurfaceFieldLimits {
            maximum_temporal_projection_terms: projection_terms - 1,
            ..exact
        })
        .is_err()
    );
    assert!(
        query(HarmonicSurfaceFieldLimits {
            maximum_quadrature_points: 0,
            ..exact
        })
        .is_err()
    );
    let shaped = HarmonicWindow::with_profile(0.0, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    let shaped_terms = mode_count(&exposure.mean, shaped).unwrap().unwrap() * (3 + 4 * 3 + 10);
    let shaped_caps = HarmonicSurfaceFieldLimits {
        maximum_temporal_projection_terms: shaped_terms,
        ..exact
    };
    exposure
        .surface_field_window(&surface(), [0.0; 3], shaped, shaped_caps, || false)
        .unwrap();
    assert!(
        exposure
            .surface_field_window(
                &surface(),
                [0.0; 3],
                shaped,
                HarmonicSurfaceFieldLimits {
                    maximum_temporal_projection_terms: shaped_terms - 1,
                    ..shaped_caps
                },
                || false
            )
            .is_err()
    );
    let broad = HarmonicWindow::with_profile(0.0, 10.0, HarmonicWindowProfile::Triangular).unwrap();
    assert!(
        exposure
            .surface_field_window(&surface(), [0.0; 3], broad, exact, || false)
            .is_err()
    );
    let invalid = [HarmonicSurfaceQuadrature {
        receiver: 1,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    assert!(
        exposure
            .surface_field_window(&invalid, [0.0; 3], window, limits(), || false)
            .is_err()
    );
    let mut callbacks = 0;
    assert!(matches!(
        exposure.surface_field_window(&surface(), [0.0; 3], window, limits(), || {
            callbacks += 1;
            callbacks == 8
        }),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(callbacks, 8);
    assert_eq!(query(exact).unwrap(), accepted);
    assert_eq!(
        exposure
            .surface_field_window(&[], [0.0; 3], window, exact, || false)
            .unwrap(),
        HarmonicSurfaceFieldMoments::default()
    );
}

#[test]
fn triangular_dark_surface_keeps_stress_and_torque_but_intervals_remain_uniform() {
    let carrier = 2.0_f64.powi(49);
    let owners = plane([carrier, carrier + 1.0, carrier + 2.0]);
    let exposure = bank(&owners, [1.0, -2.0, 1.0].map(C64::from_re), true);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let start = -0.5 * duration;
        let window =
            HarmonicWindow::with_profile(start, duration, HarmonicWindowProfile::Triangular)
                .unwrap();
        let actual = exposure
            .surface_field_window(&surface(), [0.0; 3], window, limits(), || false)
            .unwrap();
        let square = duration.powi(4) / 480.0;
        close(actual.mean_power_w, 16.0 * square, 3e-12);
        close(actual.mean_linear_momentum_rate_n[0], 4.0 * square, 3e-12);
        close(
            actual.mean_angular_momentum_rate_nm[1],
            -1.6 * square,
            3e-12,
        );
        close(
            actual.mean_angular_momentum_rate_nm[2],
            -12.0 * square,
            3e-12,
        );
        let interval = exposure
            .surface_field_interval_from_clock_dd(
                &surface(),
                [0.0; 3],
                Dd::from_f64(start),
                duration,
                limits(),
                || false,
            )
            .unwrap();
        close(interval.energy_j, 48.0 * square * duration, 3e-12);
    }
}
