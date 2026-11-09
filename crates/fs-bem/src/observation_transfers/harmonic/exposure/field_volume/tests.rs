use super::*;

fn source_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 100_000,
        maximum_receiver_frequency_pairs: 100_000,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn limits() -> HarmonicVolumeFieldLimits {
    HarmonicVolumeFieldLimits {
        maximum_quadrature_points: 128,
        maximum_field_terms: 100_000,
    }
}
fn owner(omega: f64, positions: &[[f64; 3]], fields: Vec<[C64; 4]>) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 2.0,
            sound_speed: 4.0,
        },
        sources: vec![SourceId::StructuralMode("source".into())],
        receivers: positions
            .iter()
            .enumerate()
            .map(|(i, &position_m)| ReceiverSpec {
                id: format!("physical-{i}"),
                position_m,
            })
            .collect(),
        producer_hash: ContentHash([83; 32]),
        transfers: fields,
        radiation_momentum: true,
    }
}
fn plane(frequencies: [f64; 3]) -> [PreparedHarmonicObservation; 3] {
    frequencies.map(|omega| {
        owner(
            omega,
            &[[0.0, 3.0, -0.4]],
            vec![[C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO]],
        )
    })
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
        PreparedHarmonicExposure::new(&owners[0], "source-clock", &[], &[&lines], source_limits())
            .unwrap()
    } else {
        PreparedHarmonicExposure::new(&owners[0], "source-clock", &lines, &[], source_limits())
            .unwrap()
    }
}
fn volume() -> [HarmonicVolumeQuadrature; 1] {
    [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 2.0,
    }]
}
fn close(actual: Dd, expected: f64, tolerance: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:e}"
    );
}
fn unpack(s: HarmonicVolumeFieldStorage) -> [Dd; 7] {
    [
        s.energy_j,
        s.wave_momentum.linear_momentum_ns[0],
        s.wave_momentum.linear_momentum_ns[1],
        s.wave_momentum.linear_momentum_ns[2],
        s.wave_momentum.angular_momentum_nms[0],
        s.wave_momentum.angular_momentum_nms[1],
        s.wave_momentum.angular_momentum_nms[2],
    ]
}
fn delta(s: HarmonicVolumeFieldChange) -> [Dd; 7] {
    [
        s.energy_change_j,
        s.wave_momentum.linear_momentum_ns[0],
        s.wave_momentum.linear_momentum_ns[1],
        s.wave_momentum.linear_momentum_ns[2],
        s.wave_momentum.angular_momentum_nms[0],
        s.wave_momentum.angular_momentum_nms[1],
        s.wave_momentum.angular_momentum_nms[2],
    ]
}

#[test]
fn coherent_dark_storage_and_changes_retain_proper_and_real_nulls() {
    let carrier = 2.0_f64.powi(49);
    let owners = plane([carrier, carrier + 1.0, carrier + 2.0]);
    let proper = bank(&owners, [1.0, -2.0, 1.0].map(C64::from_re), true);
    for dt in [1e-6_f64, 1e-9, 1e-12] {
        let stored = proper
            .volume_field_storage_from_clock_dd(
                &volume(),
                [0.0; 3],
                Dd::from_f64(dt),
                limits(),
                || false,
            )
            .unwrap();
        let change = proper
            .volume_field_change_from_clock_dd(&volume(), [0.0; 3], Dd::ZERO, dt, limits(), || {
                false
            })
            .unwrap();
        close(stored.energy_j, 2.0 * dt.powi(4), 3e-12);
        close(
            stored.wave_momentum.linear_momentum_ns[0],
            0.5 * dt.powi(4),
            3e-12,
        );
        close(
            stored.wave_momentum.angular_momentum_nms[2],
            -1.5 * dt.powi(4),
            3e-12,
        );
        for (a, b) in delta(change).into_iter().zip(unpack(stored)) {
            close(a, b.to_f64(), 3e-12);
        }
    }
    assert_eq!(
        proper
            .volume_field_storage_from_clock_dd(&volume(), [0.0; 3], Dd::ZERO, limits(), || false)
            .unwrap(),
        HarmonicVolumeFieldStorage::default()
    );
    let owners = plane([1.0, 2.0, 3.0]);
    let real = bank(&owners, [1.0, -2.0, 1.0].map(|x| C64::new(0.0, x)), false);
    for dt in [1e-4_f64, 1e-7, 1e-10] {
        let changed = real
            .volume_field_change_from_clock_dd(&volume(), [0.0; 3], Dd::ZERO, dt, limits(), || {
                false
            })
            .unwrap();
        // sin(t)-2sin(2t)+sin(3t) = -2t^3+O(t^5).
        close(changed.energy_change_j, 16.0 * dt.powi(6), 2e-8);
    }
}

const COEFFICIENTS: [[[f64; 2]; 3]; 3] = [
    [[0.7, 0.2], [-0.4, 0.3], [0.1, -0.6]],
    [[-0.2, 0.1], [0.3, 0.7], [0.6, -0.2]],
    [[0.4, -0.3], [0.8, 0.1], [-0.2, 0.5]],
];

fn oracle(
    owners: &[PreparedHarmonicObservation; 3],
    quadrature: &[HarmonicVolumeQuadrature],
    origin: [f64; 3],
    time: f64,
) -> [f64; 7] {
    let mut result = [0.0; 7];
    // Independent ensemble: real draw +/-1 and proper draw 1,i,-1,-i.
    for real_draw in [-1.0, 1.0] {
        for proper_draw in [
            C64::ONE,
            C64::new(0.0, 1.0),
            C64::from_re(-1.0),
            C64::new(0.0, -1.0),
        ] {
            for sample in quadrature {
                let mut value = [0.0; 4];
                for (bank, coefficients) in COEFFICIENTS.iter().enumerate() {
                    let draw = [C64::ONE, C64::from_re(real_draw), proper_draw][bank];
                    for (i, owner) in owners.iter().enumerate() {
                        let angle = owner.omega * time;
                        let coefficient = C64::new(coefficients[i][0], coefficients[i][1]);
                        let rotation = C64::new(angle.cos(), -angle.sin());
                        for c in 0..4 {
                            value[c] += (owner.transfers[sample.receiver][c]
                                * coefficient
                                * draw
                                * rotation)
                                .re;
                        }
                    }
                }
                let [p, ux, uy, uz] = value;
                let rho = owners[0].medium.density;
                let c = owners[0].medium.sound_speed;
                let weight = sample.volume_m3 / 8.0;
                result[0] +=
                    weight * 0.5 * (p * p / (rho * c * c) + rho * (ux * ux + uy * uy + uz * uz));
                let g = [ux, uy, uz].map(|u| weight * p * u / (c * c));
                let arm: [f64; 3] = core::array::from_fn(|i| {
                    owners[0].receivers[sample.receiver].position_m[i] - origin[i]
                });
                for i in 0..3 {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    result[i + 1] += g[i];
                    result[i + 4] += arm[j] * g[k] - arm[k] * g[j];
                }
            }
        }
    }
    result
}

#[test]
fn mixed_statistics_storage_and_change_match_independent_phase_ensemble() {
    let positions = [[0.2, 0.4, -0.1], [-0.3, 0.5, 0.7]];
    let owners: [_; 3] = core::array::from_fn(|i| {
        owner(
            [1.2, 1.7, 2.3][i],
            &positions,
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
        "ensemble-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        source_limits(),
    )
    .unwrap();
    let quadrature = [
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: 0.7,
        },
        HarmonicVolumeQuadrature {
            receiver: 1,
            volume_m3: 0.4,
        },
    ];
    let origin = [0.1, 0.2, -0.3];
    let first = oracle(&owners, &quadrature, origin, 0.2);
    let last = oracle(&owners, &quadrature, origin, 0.7);
    let stored = unpack(
        exposure
            .volume_field_storage_from_clock_dd(
                &quadrature,
                origin,
                Dd::from_f64(0.2),
                limits(),
                || false,
            )
            .unwrap(),
    );
    let change = delta(
        exposure
            .volume_field_change_from_clock_dd(
                &quadrature,
                origin,
                Dd::from_f64(0.2),
                0.5,
                limits(),
                || false,
            )
            .unwrap(),
    );
    for i in 0..7 {
        close(stored[i], first[i], 3e-13);
        close(change[i], last[i] - first[i], 3e-13);
    }
}

#[test]
fn analytic_changes_preserve_sub_clock_intervals_and_large_finite_storage() {
    let owners = plane([1.0, 2.0, 3.0]);
    let real = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], false);
    let dt = 1e-40;
    let start = Dd::from_pair(1.0, 1e-17);
    assert_eq!(start + Dd::from_f64(dt), start);
    let changed = real
        .volume_field_change_from_clock_dd(&volume(), [0.0; 3], start, dt, limits(), || false)
        .unwrap();
    close(changed.energy_change_j, -4.0 * (2.0_f64).sin() * dt, 3e-15);
    close(
        changed.wave_momentum.linear_momentum_ns[0],
        -(2.0_f64).sin() * dt,
        3e-15,
    );
    let epoch = Dd::from_pair(2.0_f64.powi(49), 0.001);
    let query = |t, dt| {
        delta(
            real.volume_field_change_from_clock_dd(&volume(), [0.0; 3], t, dt, limits(), || false)
                .unwrap(),
        )
    };
    let whole = query(epoch, 0.04);
    let first = query(epoch, 0.02);
    let last = query(epoch + Dd::from_f64(0.02), 0.02);
    for i in 0..7 {
        close(first[i] + last[i], whole[i].to_f64(), 5e-14);
    }
    assert!(
        (query(Dd::from_f64(epoch.to_f64()), 0.04)[0] - whole[0])
            .to_f64()
            .abs()
            > 1e-6
    );

    let mut owners = [1.0, 2.0, 3.0].map(|w| {
        owner(
            w,
            &[[0.0; 3]],
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
    let exposure = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], false);
    let q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 1e-100,
    }];
    let stored = exposure
        .volume_field_storage_from_clock_dd(&q, [0.0; 3], Dd::from_f64(0.25), limits(), || false)
        .unwrap();
    close(stored.energy_j, 2.25e300 * (0.25_f64).cos().powi(2), 4e-15);
    let changed = exposure
        .volume_field_change_from_clock_dd(
            &q,
            [0.0; 3],
            Dd::from_f64(0.25),
            1e-100,
            limits(),
            || false,
        )
        .unwrap();
    close(changed.energy_change_j, -2.25e200 * (0.5_f64).sin(), 4e-15);
    // The exact change here is approximately -2.25e-100 J, although the
    // unscaled real phase increment is below the representable DD range.
    // Refuse that range boundary rather than silently returning zero energy.
    assert!(matches!(
        exposure.volume_field_change_from_clock_dd(
            &q,
            [0.0; 3],
            Dd::ZERO,
            1e-200,
            limits(),
            || false,
        ),
        Err(TransferError::Invalid(
            "acoustic volume phase increment is unrepresentable"
        ))
    ));
}

#[test]
fn source_free_volume_and_surface_energy_and_wave_momentum_balance_converge() {
    let mut previous = f64::INFINITY;
    for n in [4usize, 8, 16, 32] {
        let length = 0.7;
        let mut positions = vec![[0.0, 0.0, 0.0], [length, 0.0, 0.0]];
        positions.extend((0..n).map(|i| [(i as f64 + 0.5) * length / n as f64, 0.0, 0.0]));
        let owners = [1.0, 2.0, 3.0].map(|w| {
            owner(
                w,
                &positions,
                positions
                    .iter()
                    .map(|p| {
                        let a = w * p[0] / 4.0;
                        let u = C64::new(a.cos(), a.sin());
                        [u * C64::from_re(8.0), u, C64::ZERO, C64::ZERO]
                    })
                    .collect(),
            )
        });
        let exposure = bank(
            &owners,
            [C64::ONE, C64::new(0.0, 0.4), C64::new(-0.25, 0.3)],
            true,
        );
        let q: Vec<_> = (0..n)
            .map(|i| HarmonicVolumeQuadrature {
                receiver: i + 2,
                volume_m3: length / n as f64,
            })
            .collect();
        let surface = [
            HarmonicSurfaceQuadrature {
                receiver: 0,
                oriented_area_m2: [-1.0, 0.0, 0.0],
            },
            HarmonicSurfaceQuadrature {
                receiver: 1,
                oriented_area_m2: [1.0, 0.0, 0.0],
            },
        ];
        let stored = exposure
            .volume_field_change_from_clock_dd(
                &q,
                [0.0; 3],
                Dd::from_f64(0.2),
                0.1,
                limits(),
                || false,
            )
            .unwrap();
        let outgoing = exposure
            .surface_field_interval_from_clock_dd(
                &surface,
                [0.0; 3],
                Dd::from_f64(0.2),
                0.1,
                HarmonicSurfaceFieldLimits {
                    maximum_quadrature_points: 2,
                    maximum_temporal_projection_terms: 10000,
                    maximum_quadrature_frequency_pairs: 1000,
                },
                || false,
            )
            .unwrap();
        let residual = (stored.energy_change_j + outgoing.energy_j).to_f64().abs();
        assert!(
            residual > 1e-12 && residual < previous * 0.26,
            "n={n}, residual={residual:e}, previous={previous:e}"
        );
        close(
            stored.wave_momentum.linear_momentum_ns[0]
                + outgoing.wave_momentum.linear_momentum_ns[0],
            (stored.energy_change_j + outgoing.energy_j).to_f64() / 4.0,
            4e-11,
        );
        previous = residual;
    }
}

#[test]
fn volume_geometry_statistics_authority_and_atomic_admission_are_retained() {
    let mut owners = plane([1.0, 2.0, 3.0]);
    let exposure = bank(&owners, [C64::ONE; 3], true);
    let q = volume();
    let query = |quadrature: &[HarmonicVolumeQuadrature], origin, caps| {
        exposure.volume_field_storage_from_clock_dd(
            quadrature,
            origin,
            Dd::from_f64(0.2),
            caps,
            || false,
        )
    };
    let exact = HarmonicVolumeFieldLimits {
        maximum_quadrature_points: 1,
        maximum_field_terms: 15,
    };
    let base = query(&q, [0.0; 3], exact).unwrap();
    assert!(matches!(
        query(
            &q,
            [0.0; 3],
            HarmonicVolumeFieldLimits {
                maximum_field_terms: 14,
                ..exact
            }
        ),
        Err(TransferError::Cap {
            requested: 15,
            cap: 14
        })
    ));
    assert!(
        query(
            &q,
            [0.0; 3],
            HarmonicVolumeFieldLimits {
                maximum_quadrature_points: 0,
                ..exact
            }
        )
        .is_err()
    );
    let split = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 1.0,
    }; 2];
    for (a, b) in unpack(query(&split, [0.0; 3], limits()).unwrap())
        .into_iter()
        .zip(unpack(base))
    {
        assert!((a - b).to_f64().abs() < 1e-28 * a.to_f64().abs().max(1.0));
    }
    let shifted = query(&q, [0.0, 1.0, 0.0], exact).unwrap();
    close(
        shifted.wave_momentum.angular_momentum_nms[2],
        (base.wave_momentum.angular_momentum_nms[2] + base.wave_momentum.linear_momentum_ns[0])
            .to_f64(),
        2e-15,
    );
    for bad_volume in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            query(
                &[HarmonicVolumeQuadrature {
                    receiver: 0,
                    volume_m3: bad_volume
                }],
                [0.0; 3],
                limits()
            )
            .is_err()
        );
    }
    assert!(
        query(
            &[HarmonicVolumeQuadrature {
                receiver: 1,
                volume_m3: 0.0
            }],
            [0.0; 3],
            limits()
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        exposure.volume_field_storage_from_clock_dd(&q, [0.0; 3], Dd::from_f64(0.2), exact, || {
            polls += 1;
            polls == 8
        }),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(polls, 8);
    assert_eq!(query(&q, [0.0; 3], exact).unwrap(), base);
    assert_eq!(
        query(&[], [0.0; 3], exact).unwrap(),
        HarmonicVolumeFieldStorage::default()
    );
    assert_eq!(
        exposure
            .volume_field_change_from_clock_dd(&q, [0.0; 3], Dd::ONE, 0.0, exact, || false)
            .unwrap(),
        HarmonicVolumeFieldChange::default()
    );
    for dt in [-0.1, f64::INFINITY, f64::NAN] {
        assert!(
            exposure
                .volume_field_change_from_clock_dd(&q, [0.0; 3], Dd::ONE, dt, exact, || false)
                .is_err()
        );
    }
    owners[2].radiation_momentum = false;
    let unauthorized = bank(&owners, [C64::ONE, C64::ZERO, C64::ZERO], true);
    assert!(
        unauthorized
            .volume_field_storage_from_clock_dd(&q, [0.0; 3], Dd::ZERO, limits(), || false)
            .is_err()
    );
    assert!(
        unauthorized
            .volume_field_change_from_clock_dd(&[], [0.0; 3], Dd::ZERO, 0.0, limits(), || false)
            .is_err()
    );
}
