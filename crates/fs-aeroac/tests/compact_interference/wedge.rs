//! Independent polar images, complete force transforms and physical wedge walls.
use super::*;
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};

fn walls(order: usize, boundaries: [IdealPlaneBoundary; 2]) -> [IdealAcousticPlane; 2] {
    let (sine, cosine) = (PI / order as f64).sin_cos();
    [
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundaries[0]).unwrap(),
        IdealAcousticPlane::new(
            [0.0; 3],
            if order == 2 {
                [1.0, 0.0, 0.0]
            } else {
                [sine, -cosine, 0.0]
            },
            boundaries[1],
        )
        .unwrap(),
    ]
}

fn point(order: usize, fraction: f64, radius: f64, z: f64) -> [f64; 3] {
    let (sine, cosine) = (fraction * PI / order as f64).sin_cos();
    [radius * cosine, radius * sine, z]
}

fn drive() -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: C64::new(2e-4, -1e-4),
        force_on_fluid_n: [C64::new(0.2, -0.3), C64::new(-0.1, 0.4), C64::new(0.6, 0.2)],
    }
}

// Direct polar group formula, independent of production's affine reflections
// and Cartesian matrix recurrence. Each source is rotated by 2 k pi/m after
// optionally replacing theta by -theta, with the same action on its force.
fn images(
    position: [f64; 3],
    input: CompactSourcePhasor,
    order: usize,
    boundaries: [IdealPlaneBoundary; 2],
) -> (Vec<[f64; 3]>, Vec<CompactSourcePhasor>) {
    let sign = |b| {
        if b == IdealPlaneBoundary::Rigid {
            1.0
        } else {
            -1.0
        }
    };
    let mut positions = Vec::new();
    let mut sources = Vec::new();
    for k in 0..order {
        let (sine, cosine) = (TAU * k as f64 / order as f64).sin_cos();
        let rotation_sign = if k % 2 == 0 {
            1.0
        } else {
            sign(boundaries[0]) * sign(boundaries[1])
        };
        for mirrored in [false, true] {
            let parity = if mirrored { -1.0 } else { 1.0 };
            let character = rotation_sign * if mirrored { sign(boundaries[0]) } else { 1.0 };
            positions.push([
                cosine * position[0] - sine * parity * position[1],
                sine * position[0] + cosine * parity * position[1],
                position[2],
            ]);
            let [fx, fy, fz] = input.force_on_fluid_n;
            sources.push(CompactSourcePhasor {
                volume_flow_m3_s: input.volume_flow_m3_s.scale(character),
                force_on_fluid_n: [
                    fx.scale(cosine) - fy.scale(sine * parity),
                    fx.scale(sine) + fy.scale(cosine * parity),
                    fz,
                ]
                .map(|v| v.scale(character)),
            });
        }
    }
    (positions, sources)
}

#[test]
fn dihedral_images_match_polar_pressure_and_full_near_field_momentum() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    for order in [2, 3, 4, 5, 6, 16] {
        let source = point(order, 0.37, 0.04, -0.02);
        let receiver = point(order, 0.61, 0.19, 0.05);
        for boundaries in [[Rigid; 2], [PressureRelease; 2], [Rigid, PressureRelease]] {
            if order % 2 != 0 && boundaries[0] != boundaries[1] {
                continue;
            }
            for frequency in [17.0, 1300.0] {
                let field = observe(1.25, 340.0, frequency, receiver, &[source])
                    .with_wedge(walls(order, boundaries), order, 2 * order)
                    .unwrap()
                    .evaluate(&[drive()])
                    .unwrap();
                let (positions, inputs) = images(source, drive(), order, boundaries);
                let pressure =
                    |x| reference_pressure(1.25, 340.0, frequency, x, &positions, &inputs);
                let expected = pressure(receiver);
                // Near release walls a complete image group can cancel most
                // digits; compare against the total independently summed scale.
                let scale: f64 = positions
                    .iter()
                    .zip(&inputs)
                    .map(|(&p, &s)| {
                        reference_pressure(1.25, 340.0, frequency, receiver, &[p], &[s]).abs()
                    })
                    .sum();
                assert!((field.pressure_pa - expected).abs() < 2e-12 * scale.max(1e-8));
                for axis in 0..3 {
                    let h = 2e-6;
                    let mut left = receiver;
                    let mut right = receiver;
                    left[axis] -= h;
                    right[axis] += h;
                    let expected = (pressure(right) - pressure(left)).scale(0.5 / h)
                        / C64::new(0.0, TAU * frequency * 1.25);
                    let actual = field.particle_velocity_m_s[axis];
                    assert!(
                        (actual - expected).abs() < 2e-6 * actual.abs().max(scale * 1e-6),
                        "m={order}, component={axis}, actual={actual:?}, expected={expected:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn both_walls_and_the_edge_hold_for_every_source_component() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    for order in [2, 3, 4, 7, 128] {
        let position = point(order, 0.37, 0.4, -0.2);
        for boundaries in [
            [Rigid; 2],
            [PressureRelease; 2],
            [Rigid, PressureRelease],
            [PressureRelease, Rigid],
        ] {
            if order % 2 != 0 && boundaries[0] != boundaries[1] {
                continue;
            }
            let planes = walls(order, boundaries);
            let n = planes[1].inward_normal();
            let receivers = [[0.8, 0.0, 0.3], [-n[1], n[0], 0.3], [0.0, 0.0, 0.8]];
            for (receiver_index, receiver) in receivers.into_iter().enumerate() {
                for component in 0..4 {
                    let mut input = CompactSourcePhasor::default();
                    if component == 0 {
                        input.volume_flow_m3_s = C64::new(1e-4, 2e-4);
                    } else {
                        input.force_on_fluid_n[component - 1] = C64::new(0.4, -0.7);
                    }
                    let direct = observe(1.25, 340.0, 317.0, receiver, &[position]);
                    let original = direct.evaluate(&[input]).unwrap();
                    let field = direct
                        .with_wedge(planes, order, 2 * order)
                        .unwrap()
                        .evaluate(&[input])
                        .unwrap();
                    for wall in 0..2 {
                        if receiver_index != wall && receiver_index != 2 {
                            continue;
                        }
                        let scale = 2.0 * order as f64;
                        match boundaries[wall] {
                            Rigid => {
                                let normal = planes[wall].inward_normal();
                                let velocity = (0..3).fold(C64::ZERO, |sum, axis| {
                                    sum + field.particle_velocity_m_s[axis].scale(normal[axis])
                                });
                                let velocity_scale: f64 =
                                    original.particle_velocity_m_s.iter().map(|v| v.abs()).sum();
                                assert!(
                                    velocity.abs() < 2e-11 * scale * velocity_scale.max(1e-10),
                                    "m={order}, wall={wall}, component={component}, residual={velocity:?}"
                                );
                            }
                            PressureRelease => assert!(
                                field.pressure_pa.abs()
                                    < 2e-11 * scale * original.pressure_pa.abs().max(1e-10)
                            ),
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn right_angle_wedge_matches_the_existing_corner_for_all_wall_parities() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let source = [0.3, 0.5, 0.7];
    let receiver = [0.8, 0.9, 1.1];
    for boundaries in [
        [Rigid; 2],
        [PressureRelease; 2],
        [Rigid, PressureRelease],
        [PressureRelease, Rigid],
    ] {
        let original = observe(1.25, 340.0, 317.0, receiver, &[source]);
        let planes = walls(2, boundaries);
        let wedge = original.clone().with_wedge(planes, 2, 4).unwrap();
        let corner = original.with_orthogonal_planes(&planes, 4).unwrap();
        let actual = wedge.evaluate(&[drive()]).unwrap();
        let expected = corner.evaluate(&[drive()]).unwrap();
        close_complex(actual.pressure_pa, expected.pressure_pa);
        for axis in 0..3 {
            close_complex(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
            );
        }
        close(
            wedge
                .evaluate_independent_sources(&[drive()])
                .unwrap()
                .mean_square_pressure_pa2,
            corner
                .evaluate_independent_sources(&[drive()])
                .unwrap()
                .mean_square_pressure_pa2,
        );
    }
}

#[test]
fn wedge_is_covariant_under_rigid_motion_and_wall_order() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let rotate = |v: [f64; 3]| [0.8 * v[0] + 0.6 * v[2], v[1], -0.6 * v[0] + 0.8 * v[2]];
    let rotate_complex = |v: [C64; 3]| {
        [
            v[0].scale(0.8) + v[2].scale(0.6),
            v[1],
            v[0].scale(-0.6) + v[2].scale(0.8),
        ]
    };
    let offset = [3.0, -2.0, 4.0];
    let move_point = |p| {
        let v = rotate(p);
        core::array::from_fn(|i| v[i] + offset[i])
    };
    for (order, boundaries) in [(3, [PressureRelease; 2]), (4, [Rigid, PressureRelease])] {
        let position = point(order, 0.37, 0.4, -0.2);
        let receiver = point(order, 0.61, 1.3, 0.8);
        let planes = walls(order, boundaries);
        let field = observe(1.25, 340.0, 317.0, receiver, &[position])
            .with_wedge(planes, order, 2 * order)
            .unwrap()
            .evaluate(&[drive()])
            .unwrap();
        let transformed = [planes[1], planes[0]].map(|p| {
            IdealAcousticPlane::new(offset, rotate(p.inward_normal()), p.boundary()).unwrap()
        });
        let input = CompactSourcePhasor {
            volume_flow_m3_s: drive().volume_flow_m3_s,
            force_on_fluid_n: rotate_complex(drive().force_on_fluid_n),
        };
        let actual = observe(
            1.25,
            340.0,
            317.0,
            move_point(receiver),
            &[move_point(position)],
        )
        .with_wedge(transformed, order, 2 * order)
        .unwrap()
        .evaluate(&[input])
        .unwrap();
        close_complex(actual.pressure_pa, field.pressure_pa);
        let expected = rotate_complex(field.particle_velocity_m_s);
        for axis in 0..3 {
            close_complex(actual.particle_velocity_m_s[axis], expected[axis]);
        }
    }
}

#[test]
fn wedge_preserves_source_coherence_cancellation_and_original_clock() {
    let order = 5;
    let planes = walls(order, [IdealPlaneBoundary::Rigid; 2]);
    let position = point(order, 0.37, 0.4, -0.2);
    let original = observe(1.25, 340.0, 317.0, [0.0, 0.0, 0.8], &[position; 3]);
    let direct = original
        .evaluate(&[
            drive(),
            CompactSourcePhasor::default(),
            CompactSourcePhasor::default(),
        ])
        .unwrap();
    let prepared = original
        .clone()
        .with_wedge(planes, order, 6 * order)
        .unwrap();
    assert_eq!(prepared.reflection_planes(), planes);
    assert_eq!(prepared.reflection_wedge_order(), Some(order));
    assert_eq!(original.reflection_wedge_order(), None);
    assert_eq!(prepared.acoustic(), original.acoustic());
    assert_eq!(prepared.source_position_m(0), Some(position));
    assert_eq!(prepared.source_delay_s(0), original.source_delay_s(0));
    assert_eq!(
        prepared.source_phase_delay_s(0),
        original.source_phase_delay_s(0)
    );
    let negative = CompactSourcePhasor {
        volume_flow_m3_s: -drive().volume_flow_m3_s,
        force_on_fluid_n: drive().force_on_fluid_n.map(|v| -v),
    };
    let inputs = [drive(), negative, CompactSourcePhasor::default()];
    assert_eq!(prepared.evaluate(&inputs).unwrap().pressure_pa, C64::ZERO);
    let independent = prepared.evaluate_independent_sources(&inputs).unwrap();
    close(
        independent.mean_square_pressure_pa2,
        (2 * order).pow(2) as f64 * direct.pressure_pa.norm_sq(),
    );
    let one = [
        drive(),
        CompactSourcePhasor::default(),
        CompactSourcePhasor::default(),
    ];
    let field = prepared.evaluate(&one).unwrap();
    let statistics = prepared
        .evaluate_coherence_modes(
            &[&one],
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: 3,
            },
        )
        .unwrap();
    close(
        statistics.mean_square_pressure_pa2,
        0.5 * field.pressure_pa.norm_sq(),
    );
    assert_eq!(statistics.scope, field.scope);
    assert!(field.scope.contains("pi/m wedge"));
    let line = prepared.evaluate_harmonic_line(&one).unwrap();
    assert_eq!(line.frequency_hz, 317.0);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    let quiet = prepared
        .evaluate(&[1e16, 1.0, -1e16].map(|q| CompactSourcePhasor::monopole(C64::from_re(q))))
        .unwrap();
    let unit = prepared
        .evaluate(&[
            CompactSourcePhasor::monopole(C64::ONE),
            CompactSourcePhasor::default(),
            CompactSourcePhasor::default(),
        ])
        .unwrap();
    close_complex(quiet.pressure_pa, unit.pressure_pa);
}

#[test]
fn molecular_wedge_retains_the_original_complex_green_function_and_momentum() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let state = relaxing_state(true);
    let order = 4;
    let boundaries = [Rigid, PressureRelease];
    let position = point(order, 0.37, 0.4, -0.2);
    let receiver = point(order, 0.61, 1.3, 0.8);
    let frequency = 900.0;
    let prepared = molecular_observer(&state, frequency, receiver, &[position])
        .with_wedge(walls(order, boundaries), order, 2 * order)
        .unwrap();
    assert_eq!(prepared.molecular_state(), Some(&state));
    let k = prepared.wavenumber_per_m();
    let (positions, sources) = images(position, drive(), order, boundaries);
    let pressure = |x: [f64; 3]| {
        positions
            .iter()
            .zip(&sources)
            .fold(C64::ZERO, |sum, (p, source)| {
                let d: [f64; 3] = core::array::from_fn(|i| x[i] - p[i]);
                let r = d.iter().map(|v| v * v).sum::<f64>().sqrt();
                let green = phase(k.re * r).scale((-k.im * r).exp() / (4.0 * PI * r));
                let force = (0..3).fold(C64::ZERO, |sum, i| {
                    sum + source.force_on_fluid_n[i].scale(d[i] / r)
                });
                sum + green
                    * (source.volume_flow_m3_s
                        * C64::new(0.0, -state.gas().density * TAU * frequency)
                        + force * C64::new(1.0 / r + k.im, -k.re))
            })
    };
    let field = prepared.evaluate(&[drive()]).unwrap();
    close_complex(field.pressure_pa, pressure(receiver));
    for axis in 0..3 {
        let mut left = receiver;
        let mut right = receiver;
        left[axis] -= 2e-6;
        right[axis] += 2e-6;
        let expected = (pressure(right) - pressure(left)).scale(1.0 / 4e-6)
            / C64::new(0.0, TAU * frequency * state.gas().density);
        assert!(
            (field.particle_velocity_m_s[axis] - expected).abs() < 3e-7 * expected.abs().max(1e-8)
        );
    }
}

#[test]
fn wedge_refuses_invalid_angle_parity_domain_budget_and_cancels() {
    use IdealPlaneBoundary::{PressureRelease, Rigid};
    let order = 3;
    let planes = walls(order, [Rigid; 2]);
    let source = point(order, 0.37, 0.4, -0.2);
    let receiver = point(order, 0.61, 1.3, 0.8);
    let original = || observe(1.25, 340.0, 317.0, receiver, &[source]);
    for invalid_order in [0, 1, 129, usize::MAX] {
        assert!(
            original()
                .with_wedge(planes, invalid_order, usize::MAX)
                .is_err()
        );
    }
    assert!(
        original()
            .with_wedge(walls(order, [Rigid, PressureRelease]), order, 6)
            .is_err()
    );
    assert!(
        original()
            .with_wedge(walls(4, [Rigid; 2]), order, 6)
            .is_err()
    );
    assert!(original().with_wedge([planes[0]; 2], order, 6).is_err());
    assert!(original().with_wedge(planes, order, 5).is_err());
    for source in [
        [0.3, 0.0, 0.2],
        [0.3, -0.1, 0.2],
        point(order, 1.1, 0.4, -0.2),
    ] {
        assert!(
            observe(1.25, 340.0, 317.0, receiver, &[source])
                .with_wedge(planes, order, 6)
                .is_err()
        );
    }
    assert!(
        observe(1.25, 340.0, 317.0, [0.3, -0.1, 0.2], &[source])
            .with_wedge(planes, order, 6)
            .is_err()
    );
    let mut total_polls = 0;
    original()
        .with_wedge_cancellable(planes, order, 6, || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, 8, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            original().with_wedge_cancellable(planes, order, 6, || {
                polls += 1;
                polls == stop
            }),
            Err(fs_aeroac::AeroacError::Cancelled { .. })
        ));
    }
    let prepared = original().with_wedge(planes, order, 6).unwrap();
    assert!(prepared.clone().with_wedge(planes, order, 6).is_err());
    assert!(prepared.clone().with_plane(planes[0], 2).is_err());
    assert!(prepared.with_orthogonal_planes(&planes[..1], 2).is_err());
    let empty = observe(1.25, 340.0, 317.0, receiver, &[])
        .with_wedge(planes, order, 0)
        .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
}
