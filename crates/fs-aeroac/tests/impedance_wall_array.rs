//! G1/G3: original Q/F coherence and wall absorption for complete source assemblies.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::full_wave::array::FullWaveCompactSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    ImpedanceWallFieldKind, PlanarImpedanceWall, PreparedImpedanceWallCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{CoherenceBudget, CompactFieldStatistics, CompactSourcePhasor};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

const HZ: f64 = 700.0;
const RHO: f64 = 1.2;
const SPEED: f64 = 343.0;
const ZERO: CompactSourcePhasor = CompactSourcePhasor {
    volume_flow_m3_s: C64::ZERO,
    force_on_fluid_n: [C64::ZERO; 3],
};
fn limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 2,
        maximum_total_evaluations: 262_144,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-9,
            maximum_evaluations: 65_536,
            maximum_depth: 36,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-9,
            absolute_velocity_tolerance_m_s_per_n: 1e-12,
            relative_tolerance: 1e-9,
            maximum_evaluations: 65_536,
            maximum_depth: 36,
        },
    }
}
fn amplitudes() -> [CompactSourcePhasor; 2] {
    [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(0.0001, -0.0003),
            force_on_fluid_n: [
                C64::new(0.2, 0.1),
                C64::new(-0.1, 0.05),
                C64::new(0.1, -0.2),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-0.0002, 0.0001),
            force_on_fluid_n: [
                C64::new(-0.04, 0.03),
                C64::new(0.06, -0.02),
                C64::new(0.08, 0.09),
            ],
        },
    ]
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= 1e-13 + tolerance * expected.abs(),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn plane(origin: [f64; 3], normal: [f64; 3], normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(RHO, SPEED).unwrap(),
        origin,
        normal,
        normalized.map(|z| z.scale(RHO * SPEED)),
    )
    .unwrap()
}
fn prepare(
    wall: PlanarImpedanceWall,
    observer: [f64; 3],
    positions: &[[f64; 3]],
    total: bool,
    controls: FullWaveCompactSourceLimits,
) -> PreparedImpedanceWallCompactSourceObserver {
    let constructor = if total {
        PreparedImpedanceWallCompactSourceObserver::new_total
    } else {
        PreparedImpedanceWallCompactSourceObserver::new_reflected
    };
    constructor(wall, HZ, observer, positions, controls, &mut || false).unwrap()
}

// Independent spherical Green function and Hessian, with physical Q/F units.
fn direct(source: [f64; 3], receiver: [f64; 3], a: CompactSourcePhasor) -> [C64; 4] {
    let d: [f64; 3] = core::array::from_fn(|i| receiver[i] - source[i]);
    let r = d.iter().map(|v| v * v).sum::<f64>().sqrt();
    let n = d.map(|v| v / r);
    let omega = TAU * HZ;
    let k = omega / SPEED;
    let (sin, cos) = (k * r).sin_cos();
    let green = C64::new(cos, sin).scale(1.0 / (4.0 * PI * r));
    let radial = green * C64::new(1.0 / r, -k);
    let derivative = green * C64::new(k * k - 2.0 / (r * r), 2.0 * k / r);
    let f = (0..3).fold(C64::ZERO, |sum, i| sum + a.force_on_fluid_n[i].scale(n[i]));
    let velocity: [C64; 3] = core::array::from_fn(|i| {
        (derivative * f.scale(n[i])
            + radial.scale(1.0 / r) * (a.force_on_fluid_n[i] - f.scale(n[i])))
            * C64::new(0.0, -1.0 / (omega * RHO))
            + radial * a.volume_flow_m3_s.scale(n[i])
    });
    [
        green * a.volume_flow_m3_s * C64::new(0.0, -omega * RHO) + radial * f,
        velocity[0],
        velocity[1],
        velocity[2],
    ]
}
fn image_field(
    wall: PlanarImpedanceWall,
    source: [f64; 3],
    receiver: [f64; 3],
    a: CompactSourcePhasor,
    sign: f64,
    total: bool,
) -> [C64; 4] {
    let n = wall.outward_normal();
    let h = (0..3)
        .map(|i| (source[i] - wall.point_m()[i]) * n[i])
        .sum::<f64>();
    let f = (0..3).fold(C64::ZERO, |sum, i| sum + a.force_on_fluid_n[i].scale(n[i]));
    let image = core::array::from_fn(|i| source[i] - 2.0 * h * n[i]);
    let reflected = direct(
        image,
        receiver,
        CompactSourcePhasor {
            volume_flow_m3_s: a.volume_flow_m3_s.scale(sign),
            force_on_fluid_n: core::array::from_fn(|i| {
                (a.force_on_fluid_n[i] - f.scale(2.0 * n[i])).scale(sign)
            }),
        },
    );
    if total {
        let incident = direct(source, receiver, a);
        core::array::from_fn(|i| incident[i] + reflected[i])
    } else {
        reflected
    }
}
fn flattened(a: [CompactSourcePhasor; 2]) -> [C64; 8] {
    core::array::from_fn(|i| {
        if i % 4 == 0 {
            a[i / 4].volume_flow_m3_s
        } else {
            a[i / 4].force_on_fluid_n[i % 4 - 1]
        }
    })
}
// Independent dense H S H* / 2, retaining all source/component covariance entries.
fn dense_moments(h: [[C64; 8]; 4], factors: &[[C64; 8]]) -> [[C64; 4]; 4] {
    let s: [[C64; 8]; 8] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            factors
                .iter()
                .fold(C64::ZERO, |sum, b| sum + b[i] * b[j].conj())
        })
    });
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut sum = C64::ZERO;
            for a in 0..8 {
                for b in 0..8 {
                    sum = sum + h[i][a] * s[a][b] * h[j][b].conj();
                }
            }
            sum.scale(0.5)
        })
    })
}
fn check_statistics(actual: &CompactFieldStatistics, expected: [[C64; 4]; 4]) {
    close(
        C64::from_re(actual.mean_square_pressure_pa2),
        expected[0][0],
        5e-11,
    );
    for i in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            expected[0][i + 1],
            5e-11,
        );
        for j in 0..3 {
            close(
                actual.velocity_covariance_m2_s2[i][j],
                expected[i + 1][j + 1],
                5e-11,
            );
        }
    }
}

#[test]
fn ideal_wall_arrays_match_rotated_analytic_q_force_images_and_dense_coherence() {
    let origin = [0.2, -0.1, 0.4];
    let normal = [0.6, 0.8, 0.0];
    let tangent = [-0.8, 0.6, 0.0];
    let positions = [
        core::array::from_fn(|i| origin[i] - 0.08 * normal[i] + 0.03 * tangent[i]),
        core::array::from_fn(|i| origin[i] - 0.06 * normal[i] - 0.02 * tangent[i]),
    ];
    let receiver = core::array::from_fn(|i| {
        origin[i] - 0.13 * normal[i] + 0.2 * tangent[i] + if i == 2 { 0.04 } else { 0.0 }
    });
    let values = amplitudes();
    let second = [
        values[1],
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(0.0003, 0.0002),
            force_on_fluid_n: [
                C64::new(0.1, -0.2),
                C64::new(0.07, 0.09),
                C64::new(-0.08, 0.03),
            ],
        },
    ];
    for (impedance, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let wall = plane(origin, normal, impedance);
        for total in [false, true] {
            let owner = prepare(
                wall,
                receiver,
                &positions,
                total,
                FullWaveCompactSourceLimits {
                    maximum_total_evaluations: 1,
                    ..limits()
                },
            );
            let mut h = [[C64::ZERO; 8]; 4];
            for column in 0..8 {
                let mut basis = ZERO;
                if column % 4 == 0 {
                    basis.volume_flow_m3_s = C64::ONE;
                } else {
                    basis.force_on_fluid_n[column % 4 - 1] = C64::ONE;
                }
                let field = image_field(wall, positions[column / 4], receiver, basis, sign, total);
                for row in 0..4 {
                    h[row][column] = field[row];
                }
            }
            let a = flattened(values);
            let expected: [C64; 4] =
                core::array::from_fn(|i| (0..8).fold(C64::ZERO, |sum, j| sum + h[i][j] * a[j]));
            let actual = owner.evaluate(&values).unwrap();
            let line = owner.evaluate_harmonic_line(&values).unwrap();
            close(actual.pressure_pa, expected[0], 5e-12);
            assert_eq!(line.frequency_hz, HZ);
            assert_eq!(line.pressure_pa, actual.pressure_pa);
            for i in 0..3 {
                close(actual.particle_velocity_m_s[i], expected[i + 1], 5e-12);
                close(
                    C64::new(
                        actual.active_intensity_w_m2[i],
                        actual.reactive_intensity_w_m2[i],
                    ),
                    expected[0].scale(0.5) * expected[i + 1].conj(),
                    1e-11,
                );
                assert_eq!(
                    line.particle_velocity_m_s[i],
                    actual.particle_velocity_m_s[i]
                );
            }
            let modes: [&[CompactSourcePhasor]; 2] = [&values, &second];
            check_statistics(
                &owner
                    .evaluate_coherence_modes(
                        &modes,
                        CoherenceBudget {
                            max_modes: 2,
                            max_source_mode_pairs: 4,
                        },
                    )
                    .unwrap(),
                dense_moments(h, &[a, flattened(second)]),
            );
            let independent = owner.evaluate_independent_sources(&values).unwrap();
            check_statistics(
                &independent,
                dense_moments(
                    h,
                    &[flattened([values[0], ZERO]), flattened([ZERO, values[1]])],
                ),
            );
            assert!(
                (independent.mean_square_pressure_pa2 - 0.5 * actual.pressure_pa.abs().powi(2))
                    .abs()
                    > 1e-4
            );
            assert_eq!(
                owner.field_kind(),
                if total {
                    ImpedanceWallFieldKind::Total
                } else {
                    ImpedanceWallFieldKind::Reflected
                }
            );
            assert_eq!(owner.source_count(), 2);
            assert_eq!(owner.source_position_m(1), Some(positions[1]));
            assert_eq!(owner.evaluations(), 0);
            assert!(
                owner
                    .source_diagnostics(0)
                    .unwrap()
                    .monopole
                    .exact_ideal_boundary
            );
            assert!(
                owner
                    .source_diagnostics(1)
                    .unwrap()
                    .force
                    .exact_ideal_boundary
            );
        }
    }
}

#[test]
fn finite_wall_array_retains_robin_absorption_and_aggregate_work_budget() {
    let wall = plane([0.0; 3], [0.0, 0.0, 1.0], Some(C64::new(0.7, 0.9)));
    let impedance = wall.impedance_pa_s_m().unwrap();
    let positions = [[-0.04, 0.02, -0.07], [0.03, -0.03, -0.05]];
    let receiver = [0.08, 0.025, 0.0];
    let values = amplitudes();
    let total = prepare(wall, receiver, &positions, true, limits());
    let reflected = prepare(wall, receiver, &positions, false, limits());
    let field = total.evaluate(&values).unwrap();
    let returning = reflected.evaluate(&values).unwrap();
    let direct0 = direct(positions[0], receiver, values[0]);
    let direct1 = direct(positions[1], receiver, values[1]);
    close(
        field.pressure_pa,
        returning.pressure_pa + direct0[0] + direct1[0],
        2e-11,
    );
    close(
        field.pressure_pa,
        impedance * field.particle_velocity_m_s[2],
        3e-8,
    );
    for i in 0..3 {
        close(
            field.particle_velocity_m_s[i],
            returning.particle_velocity_m_s[i] + direct0[i + 1] + direct1[i + 1],
            2e-11,
        );
    }
    close(
        C64::from_re(field.active_intensity_w_m2[2]),
        C64::from_re(0.5 * impedance.re * field.particle_velocity_m_s[2].abs().powi(2)),
        6e-8,
    );
    let modes: [&[CompactSourcePhasor]; 2] = [&values, &[values[1], values[0]]];
    for statistics in [
        total.evaluate_independent_sources(&values).unwrap(),
        total
            .evaluate_coherence_modes(
                &modes,
                CoherenceBudget {
                    max_modes: 2,
                    max_source_mode_pairs: 4,
                },
            )
            .unwrap(),
    ] {
        let normal_variance = statistics.velocity_covariance_m2_s2[2][2].re;
        close(
            C64::from_re(statistics.mean_square_pressure_pa2),
            C64::from_re(impedance.abs().powi(2) * normal_variance),
            6e-8,
        );
        close(
            C64::new(
                statistics.active_intensity_w_m2[2],
                statistics.reactive_intensity_w_m2[2],
            ),
            impedance.scale(normal_variance),
            6e-8,
        );
        assert!(statistics.active_intensity_w_m2[2] > 0.0);
    }
    let work = (0..2)
        .map(|i| {
            let d = total.source_diagnostics(i).unwrap();
            d.monopole.evaluations + d.force.evaluations
        })
        .sum::<usize>();
    assert_eq!(total.evaluations(), work);
    assert_eq!(reflected.evaluations(), work);
    assert!(work > 0 && work <= limits().maximum_total_evaluations);
    let bounded = FullWaveCompactSourceLimits {
        maximum_total_evaluations: work,
        ..limits()
    };
    assert_eq!(
        prepare(wall, receiver, &positions, true, bounded).evaluations(),
        work
    );
    assert!(
        PreparedImpedanceWallCompactSourceObserver::new_total(
            wall,
            HZ,
            receiver,
            &positions,
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: work - 1,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn wall_arrays_distinguish_coherent_nulls_independence_and_direct_self_fields() {
    let wall = plane([0.0; 3], [0.0, 0.0, 1.0], None);
    let position = [0.02, -0.01, -0.03];
    let a = amplitudes()[0];
    let opposite = CompactSourcePhasor {
        volume_flow_m3_s: -a.volume_flow_m3_s,
        force_on_fluid_n: a.force_on_fluid_n.map(|v| -v),
    };
    let owner = prepare(wall, position, &[position, position], false, limits());
    let silence = owner.evaluate(&[a, opposite]).unwrap();
    assert_eq!(silence.pressure_pa, C64::ZERO);
    assert_eq!(silence.particle_velocity_m_s, [C64::ZERO; 3]);
    assert!(
        owner
            .evaluate_independent_sources(&[a, opposite])
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
    assert!(
        PreparedImpedanceWallCompactSourceObserver::new_total(
            wall,
            HZ,
            position,
            &[position],
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(owner.evaluate(&[a]).is_err());
    assert!(
        owner
            .evaluate_coherence_modes(
                &[&[a, opposite]],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1
                }
            )
            .is_err()
    );
    assert!(
        owner
            .evaluate(&[a, CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0))])
            .is_err()
    );
    assert!(matches!(
        PreparedImpedanceWallCompactSourceObserver::new_total(
            wall,
            HZ,
            position,
            &[],
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let empty = prepare(wall, [0.0; 3], &[], true, limits());
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(empty.evaluations(), 0);
    assert!(
        PreparedImpedanceWallCompactSourceObserver::new_total(
            wall,
            HZ,
            [0.0, 0.0, 0.01],
            &[],
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedImpedanceWallCompactSourceObserver::new_total(
            wall,
            0.0,
            position,
            &[],
            limits(),
            &mut || false
        )
        .is_err()
    );
}
