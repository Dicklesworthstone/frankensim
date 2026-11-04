use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 64;
const DT: f64 = 1.0 / 8192.0;
const MEDIUM: Medium = Medium {
    density: 2.0,
    sound_speed: 3.0,
};
const BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 1_000_000,
};

fn measured<const S: usize>(
    segments: &[Vec<(usize, [C64; S])>],
    window: WelchWindow,
) -> WelchCrossSpectrum {
    let rows: Vec<[f64; S]> = segments
        .iter()
        .flat_map(|tones| {
            (0..N).map(move |j| {
                core::array::from_fn(|s| {
                    tones
                        .iter()
                        .map(|(bin, amplitudes)| {
                            let angle = core::f64::consts::TAU * *bin as f64 * j as f64 / N as f64;
                            amplitudes[s].re * angle.cos() + amplitudes[s].im * angle.sin()
                        })
                        .sum()
                })
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: S,
            maximum_samples_per_channel: N * segments.len(),
            maximum_transforms: S * segments.len(),
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

fn observer(
    spectrum: &WelchCrossSpectrum,
    bin: usize,
    positions: &[[f64; 3]],
    columns: &[Vec<[C64; 4]>],
) -> PreparedHarmonicObservation {
    let omega = core::f64::consts::TAU * spectrum.frequency_hz(bin).unwrap();
    let sources: Vec<_> = (0..columns.len())
        .map(|i| SourceId::PointVolumeFlow(format!("source-{i}")))
        .collect();
    let receivers: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(i, &position_m)| ReceiverSpec {
            id: format!("receiver-{i}"),
            position_m,
        })
        .collect();
    let mut samples = Vec::new();
    for column in columns {
        for field in column {
            samples.push(field[0]);
            samples.extend(
                field[1..]
                    .iter()
                    .map(|&u| C64::new(0.0, omega * MEDIUM.density) * u),
            );
        }
    }
    PreparedHarmonicObservation::from_samples(
        omega,
        MEDIUM,
        &sources,
        &receivers,
        ContentHash([11; 32]),
        samples,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn recorded_standing_waves_keep_signed_stress_and_hann_band_power() {
    let positions = [[1.0, 2.0, -1.0]];
    let columns = [
        vec![[C64::from_re(6.0), C64::ONE, C64::ZERO, C64::ZERO]],
        vec![[C64::from_re(6.0), -C64::ONE, C64::ZERO, C64::ZERO]],
    ];
    let transverse = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [0.0, 2.0, 0.0],
    }];
    for (second, expected) in [(C64::ONE, 4.0), (-C64::ONE, -4.0)] {
        for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
            let record = measured(&[vec![(3, [C64::ONE, second])]], window);
            let fields: Vec<_> = (2..=4)
                .map(|bin| observer(&record, bin, &positions, &columns))
                .collect();
            let center = fields[1]
                .momentum_surface(&transverse, [0.0; 3], || false)
                .unwrap();
            let single = center
                .evaluate_welch_bin(&record, 3, fields[1].sources(), BUDGET, || false)
                .unwrap();
            let fraction = if window == WelchWindow::PeriodicHann {
                2.0 / 3.0
            } else {
                1.0
            };
            close(single.linear_momentum_rate_n[1], expected * fraction);
            // Adding the neighboring leaked bins restores the selected tone's
            // physical power once. No extra df, source power sum or power/c.
            let lines: Vec<_> = fields
                .iter()
                .enumerate()
                .map(|(i, field)| HarmonicMomentumWelchLine { bin: i + 2, field })
                .collect();
            let band = PreparedWelchSurfaceMomentum::new(
                &record,
                &lines,
                fields[0].sources(),
                &transverse,
                [0.0; 3],
                BUDGET,
                || false,
            )
            .unwrap();
            let result = band.evaluate().unwrap();
            close(result.linear_momentum_rate_n[1], expected);
            let torque = cross(positions[0], [0.0, expected, 0.0]);
            for i in 0..3 {
                close(result.angular_momentum_rate_nm[i], torque[i]);
            }
        }
    }
    let record = measured(
        &[vec![(3, [C64::ONE, C64::ZERO])]],
        WelchWindow::Rectangular,
    );
    let field = observer(&record, 3, &positions, &columns);
    let along = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [2.0, 0.0, 0.0],
    }];
    let surface = field.momentum_surface(&along, [0.0; 3], || false).unwrap();
    let result = surface
        .evaluate_welch_bin(&record, 3, field.sources(), BUDGET, || false)
        .unwrap();
    // One plane wave: intensity=3 W/m², area=2 m², c=3 m/s.
    close(result.linear_momentum_rate_n[0], 2.0);
}

#[test]
fn broadband_full_velocity_tensor_matches_independent_time_quadrature() {
    let segments = [
        vec![
            (3, [C64::new(0.4, 0.7), C64::new(-0.3, 0.2)]),
            (9, [C64::new(0.2, -0.4), C64::new(0.7, 0.3)]),
        ],
        vec![
            (3, [C64::new(-0.2, 0.5), C64::new(0.6, -0.7)]),
            (9, [C64::new(0.5, 0.6), C64::new(-0.2, 0.1)]),
        ],
    ];
    let record = measured(&segments, WelchWindow::Rectangular);
    let positions = [[0.3, -0.5, 0.7], [-0.2, 0.4, 0.1]];
    let columns: Vec<Vec<Vec<[C64; 4]>>> = (0..2)
        .map(|f| {
            (0..2)
                .map(|s| {
                    (0..2)
                        .map(|r| {
                            core::array::from_fn(|c| {
                                C64::new(
                                    (1 + f + 2 * s + r + c) as f64
                                        * if c == 0 { 0.9 } else { 0.11 },
                                    (f as f64 - s as f64 + r as f64 - c as f64) * 0.17,
                                )
                            })
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    let fields: Vec<_> = [3, 9]
        .into_iter()
        .enumerate()
        .map(|(f, bin)| observer(&record, bin, &positions, &columns[f]))
        .collect();
    let lines = [
        HarmonicMomentumWelchLine {
            bin: 3,
            field: &fields[0],
        },
        HarmonicMomentumWelchLine {
            bin: 9,
            field: &fields[1],
        },
    ];
    let quadrature = [
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.4, -0.2, 0.5],
        },
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [-0.3, 0.7, 0.2],
        },
    ];
    let origin = [-0.4, 0.2, 0.3];
    let band = PreparedWelchSurfaceMomentum::new(
        &record,
        &lines,
        fields[0].sources(),
        &quadrature,
        origin,
        BUDGET,
        || false,
    )
    .unwrap();
    let result = band.evaluate().unwrap();
    let mut force = [0.0; 3];
    let mut torque = [0.0; 3];
    let count = 2048;
    for segment in &segments {
        for j in 0..count {
            for sample in &quadrature {
                let r = sample.receiver;
                let mut physical = [0.0; 4];
                for (f, (bin, amplitudes)) in segment.iter().enumerate() {
                    let angle =
                        core::f64::consts::TAU * *bin as f64 * (j as f64 + 0.5) / count as f64;
                    for s in 0..2 {
                        for c in 0..4 {
                            let z = columns[f][s][r][c] * amplitudes[s];
                            physical[c] += z.re * angle.cos() + z.im * angle.sin();
                        }
                    }
                }
                // Independent instantaneous scalar-potential stress followed
                // by complete-period integration, including both tones.
                let p = physical[0];
                let u = [physical[1], physical[2], physical[3]];
                let isotropic = 0.5
                    * (p * p / (MEDIUM.density * MEDIUM.sound_speed * MEDIUM.sound_speed)
                        - MEDIUM.density * u.iter().map(|x| x * x).sum::<f64>());
                let projected = u
                    .iter()
                    .zip(sample.oriented_area_m2)
                    .map(|(x, a)| x * a)
                    .sum::<f64>();
                let traction = core::array::from_fn(|i| {
                    isotropic * sample.oriented_area_m2[i] + MEDIUM.density * u[i] * projected
                });
                let arm = core::array::from_fn(|i| positions[r][i] - origin[i]);
                let moment = cross(arm, traction);
                for i in 0..3 {
                    force[i] += traction[i] / (count * segments.len()) as f64;
                    torque[i] += moment[i] / (count * segments.len()) as f64;
                }
            }
        }
    }
    for i in 0..3 {
        close(result.linear_momentum_rate_n[i], force[i]);
        close(result.angular_momentum_rate_nm[i], torque[i]);
    }
    let displacement = [0.2, -0.3, 0.5];
    let shifted = PreparedWelchSurfaceMomentum::new(
        &record,
        &lines,
        fields[0].sources(),
        &quadrature,
        core::array::from_fn(|i| origin[i] + displacement[i]),
        BUDGET,
        || false,
    )
    .unwrap()
    .evaluate()
    .unwrap();
    let correction = cross(displacement, force);
    for i in 0..3 {
        close(
            shifted.angular_momentum_rate_nm[i],
            torque[i] - correction[i],
        );
    }
}

#[test]
fn original_identity_whole_band_budget_and_late_cancellation_are_enforced() {
    let record = measured(
        &[vec![(3, [C64::ONE, C64::ONE]), (9, [C64::ONE, -C64::ONE])]],
        WelchWindow::Rectangular,
    );
    let columns = [vec![[C64::ONE; 4]], vec![[C64::new(0.3, 0.2); 4]]];
    let positions = [[0.0; 3]];
    let first = observer(&record, 3, &positions, &columns);
    let second = observer(&record, 9, &positions, &columns);
    let lines = [
        HarmonicMomentumWelchLine {
            bin: 3,
            field: &first,
        },
        HarmonicMomentumWelchLine {
            bin: 9,
            field: &second,
        },
    ];
    let quadrature = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [1.0, 0.0, 0.0],
    }];
    let band = PreparedWelchSurfaceMomentum::new(
        &record,
        &lines,
        first.sources(),
        &quadrature,
        [0.0; 3],
        BUDGET,
        || false,
    )
    .unwrap();
    assert_eq!(band.work_products(), 2 * (2 * 2 + 4 * 2 + 7));
    let exact = HarmonicWorkBudget {
        maximum_products: band.work_products(),
    };
    assert!(
        PreparedWelchSurfaceMomentum::new(
            &record,
            &lines,
            first.sources(),
            &quadrature,
            [0.0; 3],
            exact,
            || false
        )
        .is_ok()
    );
    assert!(
        PreparedWelchSurfaceMomentum::new(
            &record,
            &lines,
            first.sources(),
            &quadrature,
            [0.0; 3],
            HarmonicWorkBudget {
                maximum_products: exact.maximum_products - 1
            },
            || false
        )
        .is_err()
    );
    let expected = band.evaluate().unwrap();
    let mut calls = 0;
    assert!(matches!(
        band.evaluate_cancellable(|| {
            calls += 1;
            calls == 7
        }),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(band.evaluate().unwrap(), expected);
    let reversed = [lines[1], lines[0]];
    let repeated = [lines[0], lines[0]];
    for invalid in [&reversed[..], &repeated[..], &[]] {
        assert!(
            PreparedWelchSurfaceMomentum::new(
                &record,
                invalid,
                first.sources(),
                &quadrature,
                [0.0; 3],
                BUDGET,
                || false
            )
            .is_err()
        );
    }
    let mut changed = observer(&record, 9, &positions, &columns);
    changed.producer_hash = ContentHash([12; 32]);
    let bad = [
        lines[0],
        HarmonicMomentumWelchLine {
            bin: 9,
            field: &changed,
        },
    ];
    assert!(
        PreparedWelchSurfaceMomentum::new(
            &record,
            &bad,
            first.sources(),
            &quadrature,
            [0.0; 3],
            BUDGET,
            || false
        )
        .is_err()
    );
    let mut changed = observer(&record, 9, &positions, &columns);
    changed.radiation_momentum = false;
    let bad = [
        lines[0],
        HarmonicMomentumWelchLine {
            bin: 9,
            field: &changed,
        },
    ];
    assert!(
        PreparedWelchSurfaceMomentum::new(
            &record,
            &bad,
            first.sources(),
            &quadrature,
            [0.0; 3],
            BUDGET,
            || false
        )
        .is_err()
    );
    let reversed_sources = [first.sources()[1].clone(), first.sources()[0].clone()];
    assert!(
        PreparedWelchSurfaceMomentum::new(
            &record,
            &lines,
            &reversed_sources,
            &quadrature,
            [0.0; 3],
            BUDGET,
            || false
        )
        .is_err()
    );
    let surface = first
        .momentum_surface(&quadrature, [0.0; 3], || false)
        .unwrap();
    for bin in [0, 9, N / 2] {
        assert!(
            surface
                .evaluate_welch_bin(&record, bin, first.sources(), BUDGET, || false)
                .is_err()
        );
    }
}

#[test]
fn source_cancellation_precedes_stress_and_preserves_a_weak_residual() {
    let record = measured(&[vec![(3, [C64::ONE; 3])]], WelchWindow::Rectangular);
    let right = [C64::from_re(6.0), C64::ONE, C64::ZERO, C64::ZERO];
    let columns: Vec<_> = [1e16, 1.0, -1e16]
        .into_iter()
        .map(|s| vec![right.map(|z| z.scale(s))])
        .collect();
    let field = observer(&record, 3, &[[0.0; 3]], &columns);
    let quadrature = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [2.0, 0.0, 0.0],
    }];
    let surface = field
        .momentum_surface(&quadrature, [0.0; 3], || false)
        .unwrap();
    let result = surface
        .evaluate_welch_bin(&record, 3, field.sources(), BUDGET, || false)
        .unwrap();
    close(result.linear_momentum_rate_n[0], 2.0);
    assert_eq!(result.angular_momentum_rate_nm, [0.0; 3]);
}
