use super::*;
use crate::helmholtz::Medium;

fn medium() -> Medium {
    Medium {
        density: 2.0,
        sound_speed: 4.0,
    }
}

fn limits() -> HarmonicMomentumLimits {
    HarmonicMomentumLimits {
        maximum_quadrature_points: 128,
        maximum_quadrature_frequency_pairs: 100_000,
        maximum_generated_terms: 100_000,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn exposure_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 100_000,
        maximum_receiver_frequency_pairs: 100_000,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn observer(
    omega: f64,
    positions: &[[f64; 3]],
    columns: &[Vec<[C64; 4]>],
) -> PreparedHarmonicObservation {
    let medium = medium();
    let sources: Vec<_> = (0..columns.len())
        .map(|s| SourceId::StructuralMode(format!("source-{s}")))
        .collect();
    let receivers: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(r, &position_m)| ReceiverSpec {
            id: format!("receiver-{r}"),
            position_m,
        })
        .collect();
    let mut samples = Vec::new();
    for column in columns {
        assert_eq!(column.len(), positions.len());
        for field in column {
            samples.push(field[0]);
            samples.extend(
                field[1..]
                    .iter()
                    .map(|&u| C64::new(0.0, omega * medium.density) * u),
            );
        }
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        medium,
        &sources,
        &receivers,
        ContentHash([23; 32]),
        samples,
    )
    .unwrap()
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

fn stored(x: HarmonicWaveMomentumStorage) -> [Dd; 6] {
    [
        x.linear_momentum_ns[0],
        x.linear_momentum_ns[1],
        x.linear_momentum_ns[2],
        x.angular_momentum_nms[0],
        x.angular_momentum_nms[1],
        x.angular_momentum_nms[2],
    ]
}

fn transfer(x: HarmonicWaveMomentumTransfer) -> [Dd; 6] {
    [
        x.linear_momentum_ns[0],
        x.linear_momentum_ns[1],
        x.linear_momentum_ns[2],
        x.angular_momentum_nms[0],
        x.angular_momentum_nms[1],
        x.angular_momentum_nms[2],
    ]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance={tolerance:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn fixture(omega: f64) -> PreparedHarmonicObservation {
    observer(
        omega,
        &[[0.2, 0.5, -0.3], [-0.4, 0.1, 0.7]],
        &[
            vec![
                [
                    C64::new(1.2, 0.3),
                    C64::new(0.1, -0.2),
                    C64::new(-0.3, 0.04),
                    C64::new(0.02, 0.1),
                ],
                [
                    C64::new(-0.3, 0.5),
                    C64::new(0.04, 0.2),
                    C64::new(0.1, 0.04),
                    C64::new(-0.2, 0.1),
                ],
            ],
            vec![
                [
                    C64::new(-0.2, 0.7),
                    C64::new(0.2, 0.04),
                    C64::new(0.1, -0.2),
                    C64::new(-0.1, 0.3),
                ],
                [
                    C64::new(0.3, -0.4),
                    C64::new(-0.2, 0.1),
                    C64::new(0.04, 0.2),
                    C64::new(0.1, 0.2),
                ],
            ],
        ],
    )
}

fn quadrature() -> [HarmonicVolumeQuadrature; 2] {
    [
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: 0.7,
        },
        HarmonicVolumeQuadrature {
            receiver: 1,
            volume_m3: 1.3,
        },
    ]
}

// Independent instantaneous real-field summation. No product spectrum,
// exposure contraction, DD temporal kernel or production volume helper.
fn real_bank(lines: &[HarmonicExposureLine<'_>], r: usize, t: f64, factor: C64) -> [f64; 4] {
    let mut values = [0.0; 4];
    for line in lines {
        let (sin, cos) = (line.observation.omega * t).sin_cos();
        for (s, coefficient) in line.coefficients.iter().enumerate() {
            let field = line.observation.transfers[s * line.observation.receivers.len() + r];
            for i in 0..4 {
                let a = field[i] * *coefficient * factor;
                values[i] += a.re * cos + a.im * sin;
            }
        }
    }
    values
}

fn oracle(
    mean: &[HarmonicExposureLine<'_>],
    factors: &[(
        &[HarmonicExposureLine<'_>],
        HarmonicExposureFactorStatistics,
    )],
    reference: &PreparedHarmonicObservation,
    q: &[HarmonicVolumeQuadrature],
    origin: [f64; 3],
    time: f64,
) -> [f64; 6] {
    let states: usize = factors
        .iter()
        .map(|(_, law)| {
            if *law == HarmonicExposureFactorStatistics::Real {
                2
            } else {
                4
            }
        })
        .product();
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let mut result = [0.0; 6];
    for state in 0..states {
        for sample in q {
            let mut field = real_bank(mean, sample.receiver, time, C64::ONE);
            let mut code = state;
            for (bank, law) in factors {
                let radix = if *law == HarmonicExposureFactorStatistics::Real {
                    2
                } else {
                    4
                };
                let add = real_bank(bank, sample.receiver, time, phases[code % radix]);
                code /= radix;
                for i in 0..4 {
                    field[i] += add[i];
                }
            }
            let c = reference.medium.sound_speed;
            let linear =
                core::array::from_fn(|i| sample.volume_m3 * field[0] * field[i + 1] / (c * c));
            let arm = core::array::from_fn(|i| {
                reference.receivers[sample.receiver].position_m[i] - origin[i]
            });
            let angular = cross(arm, linear);
            for i in 0..3 {
                result[i] += linear[i] / states as f64;
                result[i + 3] += angular[i] / states as f64;
            }
        }
    }
    result
}

/// G1: independently enumerate the declared real/proper ensemble and sum the
/// original instantaneous p*u/c² and its moment at each authored volume point.
#[test]
fn storage_matches_independent_mixed_real_proper_fields_and_endpoint_changes() {
    let owners = [fixture(7.0), fixture(11.0)];
    let a = [C64::new(0.7, -0.2), C64::new(-0.1, 0.3)];
    let b = [C64::new(-0.3, 0.4), C64::new(0.2, 0.1)];
    let d = [C64::new(0.1, 0.5), C64::new(-0.4, 0.2)];
    let mean = [line(&owners[0], &a), line(&owners[1], &b)];
    let real = [line(&owners[0], &d), line(&owners[1], &a)];
    let proper = [line(&owners[0], &b), line(&owners[1], &d)];
    let laws = [
        HarmonicExposureFactorStatistics::Real,
        HarmonicExposureFactorStatistics::ProperComplex,
    ];
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "original-clock",
        &mean,
        &[&real, &proper],
        &laws,
        exposure_limits(),
    )
    .unwrap();
    let q = quadrature();
    let origin = [0.1, -0.2, 0.3];
    let volume = PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, origin, limits()).unwrap();
    assert_eq!(volume.factor_statistics(), &laws);
    assert_eq!(volume.source_clock(), "original-clock");
    assert_eq!(volume.producer_hash(), owners[0].producer_hash());
    let factors = [(&real[..], laws[0]), (&proper[..], laws[1])];
    for time in [0.0, 0.17, 0.83] {
        let expected = oracle(&mean, &factors, &owners[0], &q, origin, time);
        for (actual, expected) in stored(volume.stored_momentum_dd(time).unwrap())
            .into_iter()
            .zip(expected)
        {
            close(actual.to_f64(), expected, 3e-14);
        }
    }
    let (start, duration) = (0.17, 0.29);
    let before = oracle(&mean, &factors, &owners[0], &q, origin, start);
    let after = oracle(&mean, &factors, &owners[0], &q, origin, start + duration);
    let actual = transfer(volume.momentum_change_dd(start, duration).unwrap());
    for i in 0..6 {
        close(actual[i].to_f64(), after[i] - before[i], 5e-14);
    }
}

/// G1: a source-free two-frequency plane wave crosses an actual rectangular
/// control volume. Six outward face weights are exact for its surface flux;
/// midpoint volume cubature is deliberately approximate. Storage change plus
/// transport converges at second order under 4/8/16/32-cell refinement.
#[test]
fn plane_wave_storage_and_outward_transport_balance_under_volume_refinement() {
    let bounds = [[-0.2, 0.7], [0.3, 0.9], [-0.5, 0.4]];
    let lengths: [f64; 3] = core::array::from_fn(|i| bounds[i][1] - bounds[i][0]);
    let center: [f64; 3] = core::array::from_fn(|i| 0.5 * (bounds[i][0] + bounds[i][1]));
    let origin = [-0.2, 0.1, -0.4];
    let coefficients = [[C64::new(0.7, 0.2)], [C64::new(-0.3, 0.4)]];
    let frequencies = [5.0, 11.0];
    let mut errors = Vec::new();
    let mut original_flux = None;
    for cells in [4, 8, 16, 32] {
        let mut positions = Vec::new();
        let mut surface_q = Vec::new();
        for axis in 0..3 {
            let area = lengths[(axis + 1) % 3] * lengths[(axis + 2) % 3];
            for side in 0..2 {
                let mut position = center;
                position[axis] = bounds[axis][side];
                let mut normal_area = [0.0; 3];
                normal_area[axis] = if side == 0 { -area } else { area };
                surface_q.push(HarmonicSurfaceQuadrature {
                    receiver: positions.len(),
                    oriented_area_m2: normal_area,
                });
                positions.push(position);
            }
        }
        let dx = lengths[0] / cells as f64;
        let mut volume_q = Vec::new();
        for cell in 0..cells {
            let mut position = center;
            position[0] = bounds[0][0] + dx * (cell as f64 + 0.5);
            volume_q.push(HarmonicVolumeQuadrature {
                receiver: positions.len(),
                volume_m3: dx * lengths[1] * lengths[2],
            });
            positions.push(position);
        }
        let owners: Vec<_> = frequencies
            .iter()
            .map(|&omega| {
                let fields = positions
                    .iter()
                    .map(|p| {
                        let (sin, cos) = (omega * p[0] / medium().sound_speed).sin_cos();
                        let pressure = C64::new(cos, sin);
                        [
                            pressure,
                            pressure.scale(1.0 / (medium().density * medium().sound_speed)),
                            C64::ZERO,
                            C64::ZERO,
                        ]
                    })
                    .collect();
                observer(omega, &positions, &[fields])
            })
            .collect();
        let lines = [
            line(&owners[0], &coefficients[0]),
            line(&owners[1], &coefficients[1]),
        ];
        let exposure =
            PreparedHarmonicExposure::new(&owners[0], "plane-wave", &lines, &[], exposure_limits())
                .unwrap();
        let surface =
            PreparedHarmonicSurfaceWaveMomentum::new(&exposure, &surface_q, origin, limits())
                .unwrap();
        let volume =
            PreparedHarmonicVolumeWaveMomentum::new(&exposure, &volume_q, origin, limits())
                .unwrap();
        let outward = surface.transported_momentum_dd(0.13, 0.17).unwrap();
        let change = volume.momentum_change_dd(0.13, 0.17).unwrap();
        if let Some(first) = original_flux {
            assert_eq!(outward, first);
        } else {
            original_flux = Some(outward);
        }
        let f = transfer(outward);
        let d = transfer(change);
        assert!(f[0].to_f64().abs() > 3e-4);
        let error = (d[0] + f[0]).to_f64();
        assert!(
            error > 0.0,
            "the volume approximation must not fake exact cancellation"
        );
        close((d[1] + f[1]).to_f64(), 0.0, 1e-16);
        close((d[2] + f[2]).to_f64(), 0.0, 1e-16);
        close((d[3] + f[3]).to_f64(), 0.0, 1e-16);
        close(
            (d[4] + f[4]).to_f64(),
            (center[2] - origin[2]) * error,
            2e-16,
        );
        close(
            (d[5] + f[5]).to_f64(),
            -(center[1] - origin[1]) * error,
            2e-16,
        );
        errors.push(error);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.8..4.5).contains(&ratio), "cubature errors {errors:?}");
    }
    assert!(errors[3] < 1.3e-7, "finest cubature error {errors:?}");
}

/// G3: pressure is scalar; velocity, geometry, linear momentum and angular
/// momentum rotate together. Origin shifts retain the explicit orbital term.
#[test]
fn stored_momentum_rotates_and_changes_origin_consistently() {
    let owner = fixture(7.0);
    let coefficients = [C64::new(0.7, -0.2), C64::new(-0.1, 0.4)];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &lines, &[&lines], exposure_limits())
            .unwrap();
    let q = quadrature();
    let origin = [0.1, -0.3, 0.2];
    let volume = PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, origin, limits()).unwrap();
    let v = volume.stored_momentum_dd(0.19).unwrap();
    let shift = [0.3, -0.5, 0.1];
    let moved_origin = core::array::from_fn(|i| origin[i] + shift[i]);
    let moved = PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, moved_origin, limits())
        .unwrap()
        .stored_momentum_dd(0.19)
        .unwrap();
    let orbital = cross(shift, v.linear_momentum_ns.map(Dd::to_f64));
    for i in 0..3 {
        close(
            (moved.linear_momentum_ns[i] - v.linear_momentum_ns[i]).to_f64(),
            0.0,
            2e-15,
        );
        close(
            (moved.angular_momentum_nms[i] - v.angular_momentum_nms[i]).to_f64(),
            -orbital[i],
            2e-15,
        );
    }
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let positions: Vec<_> = owner
        .receivers
        .iter()
        .map(|r| rotate(r.position_m))
        .collect();
    let columns: Vec<_> = (0..owner.sources.len())
        .map(|s| {
            (0..owner.receivers.len())
                .map(|r| {
                    let f = owner.transfers[s * owner.receivers.len() + r];
                    [f[0], f[3], f[1], f[2]]
                })
                .collect()
        })
        .collect();
    let rotated = observer(owner.omega, &positions, &columns);
    let rotated_lines = [line(&rotated, &coefficients)];
    let rotated_exposure = PreparedHarmonicExposure::new(
        &rotated,
        "clock",
        &rotated_lines,
        &[&rotated_lines],
        exposure_limits(),
    )
    .unwrap();
    let rotated_volume =
        PreparedHarmonicVolumeWaveMomentum::new(&rotated_exposure, &q, rotate(origin), limits())
            .unwrap();
    let actual = rotated_volume.stored_momentum_dd(0.19).unwrap();
    for (a, b) in actual
        .linear_momentum_ns
        .into_iter()
        .zip(rotate(v.linear_momentum_ns.map(Dd::to_f64)))
    {
        close(a.to_f64(), b, 3e-15);
    }
    for (a, b) in actual
        .angular_momentum_nms
        .into_iter()
        .zip(rotate(v.angular_momentum_nms.map(Dd::to_f64)))
    {
        close(a.to_f64(), b, 3e-15);
    }
}

/// G0/G3: original exposure banks can carry enormous independent DC storage
/// alongside a weak varying mean. Storage increments retain that variation
/// without subtracting rounded DC endpoints or collapsing the DD clock.
#[test]
fn weak_storage_changes_survive_large_dc_and_dd_clock_subdivision() {
    let owner = observer(
        1.0,
        &[[0.0, 2.0, 0.0]],
        &[vec![[C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO]]],
    );
    let small = [C64::ONE];
    let huge = [C64::from_re(1e100)];
    let mean = [line(&owner, &small)];
    let factor = [line(&owner, &huge)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &mean, &[&factor], exposure_limits())
            .unwrap();
    let only_mean =
        PreparedHarmonicExposure::new(&owner, "clock", &mean, &[], exposure_limits()).unwrap();
    let q = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 1.0,
    }];
    let volume =
        PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, [0.0; 3], limits()).unwrap();
    let weak = PreparedHarmonicVolumeWaveMomentum::new(&only_mean, &q, [0.0; 3], limits()).unwrap();
    let start = Dd::from_pair(2.0_f64.powi(54), 0.25);
    let duration = 0.125;
    assert_eq!(start.hi + duration, start.hi);
    let before = volume.stored_momentum_from_clock_dd(start).unwrap();
    let after = volume
        .stored_momentum_from_clock_dd(start + Dd::from_f64(duration))
        .unwrap();
    assert_eq!(
        before.linear_momentum_ns[0].to_f64(),
        after.linear_momentum_ns[0].to_f64()
    );
    let whole = volume
        .momentum_change_from_clock_dd(start, duration)
        .unwrap();
    assert_ne!(whole.linear_momentum_ns[0], Dd::ZERO);
    assert_eq!(
        whole,
        weak.momentum_change_from_clock_dd(start, duration).unwrap()
    );
    let first = transfer(volume.momentum_change_from_clock_dd(start, 0.0625).unwrap());
    let last = transfer(
        volume
            .momentum_change_from_clock_dd(start + Dd::from_f64(0.0625), 0.0625)
            .unwrap(),
    );
    for (i, value) in transfer(whole).into_iter().enumerate() {
        close((value - first[i] - last[i]).to_f64(), 0.0, 2e-16);
    }
    close(
        (whole.angular_momentum_nms[2] + Dd::from_f64(2.0) * whole.linear_momentum_ns[0]).to_f64(),
        0.0,
        2e-16,
    );
}

/// G0/G4/G5: positive-volume admission, original medium authority, exact zero
/// entries and bounded cancellation are checked before publication.
#[test]
fn volume_admission_zero_weights_and_cancellation_preserve_replay() {
    let owner = fixture(7.0);
    let coefficients = [C64::ONE, C64::ZERO];
    let lines = [line(&owner, &coefficients)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &lines, &[], exposure_limits()).unwrap();
    let q = quadrature();
    let volume =
        PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, [0.0; 3], limits()).unwrap();
    for budget in [
        HarmonicMomentumLimits {
            maximum_quadrature_points: q.len() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_quadrature_frequency_pairs: volume.quadrature_frequency_pairs() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_generated_terms: volume.generated_terms() - 1,
            ..limits()
        },
        HarmonicMomentumLimits {
            maximum_workspace_bytes: volume.workspace_bytes() - 1,
            ..limits()
        },
    ] {
        assert!(matches!(
            PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, [0.0; 3], budget),
            Err(TransferError::Cap { .. })
        ));
    }
    for bad in [
        HarmonicVolumeQuadrature {
            receiver: usize::MAX,
            volume_m3: 0.0,
        },
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: -1.0,
        },
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: f64::NAN,
        },
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: f64::INFINITY,
        },
    ] {
        assert!(
            PreparedHarmonicVolumeWaveMomentum::new(&exposure, &[bad], [0.0; 3], limits()).is_err()
        );
    }
    assert!(
        PreparedHarmonicVolumeWaveMomentum::new(&exposure, &q, [f64::NAN, 0.0, 0.0], limits())
            .is_err()
    );
    let no_stress = fixture(11.0).without_radiation_momentum();
    let zero_coefficients = [C64::ZERO; 2];
    let bad_lines = [line(&no_stress, &zero_coefficients)];
    let unsupported =
        PreparedHarmonicExposure::new(&owner, "clock", &bad_lines, &[], exposure_limits()).unwrap();
    assert!(
        PreparedHarmonicVolumeWaveMomentum::new(&unsupported, &[], [0.0; 3], limits()).is_err()
    );
    let dark = PreparedHarmonicExposure::new(&owner, "dark", &[], &[], exposure_limits()).unwrap();
    let dark_volume =
        PreparedHarmonicVolumeWaveMomentum::new(&dark, &q, [0.0; 3], limits()).unwrap();
    assert_eq!(
        dark_volume.stored_momentum_dd(0.1).unwrap(),
        HarmonicWaveMomentumStorage::default()
    );
    assert_eq!(
        volume.momentum_change_dd(0.1, 0.0).unwrap(),
        HarmonicWaveMomentumTransfer::default()
    );
    assert!(volume.stored_momentum_dd(f64::NAN).is_err());
    assert!(volume.momentum_change_dd(0.0, -0.1).is_err());
    let distant = observer(
        1.0,
        &[[f64::MAX, 0.0, 0.0]],
        &[vec![[C64::from_re(8.0), C64::ONE, C64::ZERO, C64::ZERO]]],
    );
    let one = [C64::ONE];
    let distant_lines = [line(&distant, &one)];
    let distant_exposure =
        PreparedHarmonicExposure::new(&distant, "clock", &distant_lines, &[], exposure_limits())
            .unwrap();
    let zero = [HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 0.0,
    }];
    let zero_volume = PreparedHarmonicVolumeWaveMomentum::new(
        &distant_exposure,
        &zero,
        [-f64::MAX, 0.0, 0.0],
        limits(),
    )
    .unwrap();
    assert_eq!(zero_volume.quadrature_frequency_pairs(), 1);
    assert_eq!(
        zero_volume.stored_momentum_dd(0.2).unwrap(),
        HarmonicWaveMomentumStorage::default()
    );
    assert_eq!(
        zero_volume.momentum_change_dd(0.2, 0.1).unwrap(),
        HarmonicWaveMomentumTransfer::default()
    );
    let mut polls = 0;
    PreparedHarmonicVolumeWaveMomentum::new_cancellable(&exposure, &q, [0.0; 3], limits(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            PreparedHarmonicVolumeWaveMomentum::new_cancellable(
                &exposure,
                &q,
                [0.0; 3],
                limits(),
                || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(calls, stop);
    }
    let start = Dd::from_f64(0.13);
    let expected = volume.momentum_change_from_clock_dd(start, 0.29).unwrap();
    let mut change_polls = 0;
    volume
        .momentum_change_from_clock_dd_cancellable(start, 0.29, || {
            change_polls += 1;
            false
        })
        .unwrap();
    let mut endpoint_polls = 0;
    let endpoint = volume
        .stored_momentum_from_clock_dd_cancellable(start, || {
            endpoint_polls += 1;
            false
        })
        .unwrap();
    for (endpoint_query, count) in [(false, change_polls), (true, endpoint_polls)] {
        for stop in [1, count / 2, count] {
            let mut calls = 0;
            let cancelled = if endpoint_query {
                matches!(
                    volume.stored_momentum_from_clock_dd_cancellable(start, || {
                        calls += 1;
                        calls == stop
                    }),
                    Err(TransferError::Cancelled)
                )
            } else {
                matches!(
                    volume.momentum_change_from_clock_dd_cancellable(start, 0.29, || {
                        calls += 1;
                        calls == stop
                    }),
                    Err(TransferError::Cancelled)
                )
            };
            assert!(cancelled);
            assert_eq!(calls, stop);
            assert_eq!(
                volume.momentum_change_from_clock_dd(start, 0.29).unwrap(),
                expected
            );
            assert_eq!(
                volume.stored_momentum_from_clock_dd(start).unwrap(),
                endpoint
            );
        }
    }
}
