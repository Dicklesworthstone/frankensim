//! G1/G3/G4: physical tensor images, boundary identities and prepared budgets.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::{CoherenceBudget, IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn owner(frequency: f64, receiver: [f64; 3], positions: &[[f64; 3]]) -> PreparedQuadrupoleObserver {
    PreparedQuadrupoleObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        frequency,
        receiver,
        positions,
        positions.len(),
    )
    .unwrap()
}

fn source() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.4, 0.2),
            C64::new(-0.3, 0.1),
            C64::new(0.2, -0.1),
            C64::new(0.5, -0.2),
            C64::new(-0.1, 0.4),
            C64::new(0.3, 0.1),
        ],
    }
}

fn close(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * scale.max(1e-20),
        "actual={actual:?}, expected={expected:?}, scale={scale}"
    );
}

#[test]
fn g1_oblique_wall_matches_independently_reflected_tensor_including_shear() {
    let normal = [0.6, 0.8, 0.0];
    let point = [0.2, -0.1, 0.3];
    let position = [0.4, 0.5, 0.1];
    let receiver = [0.7, 0.3, 0.6];
    let q = source();
    let tensor = [
        [
            q.stress_moment_n_m[0],
            q.stress_moment_n_m[3],
            q.stress_moment_n_m[4],
        ],
        [
            q.stress_moment_n_m[3],
            q.stress_moment_n_m[1],
            q.stress_moment_n_m[5],
        ],
        [
            q.stress_moment_n_m[4],
            q.stress_moment_n_m[5],
            q.stress_moment_n_m[2],
        ],
    ];
    let h: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| (if i == j { 1.0 } else { 0.0 }) - 2.0 * normal[i] * normal[j])
    });
    let distance: f64 = (0..3).map(|i| normal[i] * (position[i] - point[i])).sum();
    let image = core::array::from_fn(|i| position[i] - 2.0 * distance * normal[i]);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let sign = if boundary == IdealPlaneBoundary::Rigid {
            1.0
        } else {
            -1.0
        };
        let transformed: [[C64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let mut sum = C64::ZERO;
                for a in 0..3 {
                    for b in 0..3 {
                        sum = sum + tensor[a][b].scale(sign * h[i][a] * h[j][b]);
                    }
                }
                sum
            })
        });
        let image_q = QuadrupolePhasor {
            stress_moment_n_m: [
                transformed[0][0],
                transformed[1][1],
                transformed[2][2],
                transformed[0][1],
                transformed[0][2],
                transformed[1][2],
            ],
        };
        let wall = IdealAcousticPlane::new(point, normal, boundary).unwrap();
        for frequency in [0.5, 127.0, 2500.0] {
            let direct = owner(frequency, receiver, &[position]);
            let reflected = direct.clone().with_orthogonal_planes(&[wall], 2).unwrap();
            let actual = reflected.evaluate(&[q]).unwrap();
            let expected = owner(frequency, receiver, &[position, image])
                .evaluate(&[q, image_q])
                .unwrap();
            close(
                actual.pressure_pa,
                expected.pressure_pa,
                expected.pressure_pa.abs(),
            );
            let velocity_scale = expected
                .particle_velocity_m_s
                .iter()
                .fold(0.0_f64, |s, u| s.max(u.abs()));
            for i in 0..3 {
                close(
                    actual.particle_velocity_m_s[i],
                    expected.particle_velocity_m_s[i],
                    velocity_scale,
                );
            }
            assert_eq!(reflected.source_count(), 1);
            assert_eq!(reflected.source_position_m(0), direct.source_position_m(0));
            assert_eq!(reflected.source_delay_s(0), direct.source_delay_s(0));
            let line = reflected.evaluate_harmonic_line(&[q]).unwrap();
            assert_eq!(line.pressure_pa, actual.pressure_pa);
            let stats = reflected
                .evaluate_coherence_modes(
                    &[&[q]],
                    CoherenceBudget {
                        max_modes: 1,
                        max_source_mode_pairs: 1,
                    },
                )
                .unwrap();
            assert_eq!(stats.scope, reflected.scope());
            close(
                C64::from_re(stats.mean_square_pressure_pa2),
                C64::from_re(0.5 * actual.pressure_pa.norm_sq()),
                actual.pressure_pa.norm_sq(),
            );
        }
    }
}

#[test]
fn g1_three_wall_mixed_corner_enforces_each_boundary_for_every_tensor_component() {
    let boundaries = [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
        IdealPlaneBoundary::Rigid,
    ];
    let planes: [_; 3] = core::array::from_fn(|axis| {
        let mut normal = [0.0; 3];
        normal[axis] = 1.0;
        IdealAcousticPlane::new([0.0; 3], normal, boundaries[axis]).unwrap()
    });
    for face in 0..3 {
        let mut receiver = [0.5, 0.4, 0.3];
        receiver[face] = 0.0;
        let direct = owner(73.0, receiver, &[[0.2, 0.1, 0.4]]);
        let bounded = direct.clone().with_orthogonal_planes(&planes, 8).unwrap();
        for component in 0..6 {
            let mut q = QuadrupolePhasor::default();
            q.stress_moment_n_m[component] = C64::ONE;
            let reference = direct.evaluate(&[q]).unwrap();
            let field = bounded.evaluate(&[q]).unwrap();
            if boundaries[face] == IdealPlaneBoundary::Rigid {
                close(
                    field.particle_velocity_m_s[face],
                    C64::ZERO,
                    reference
                        .particle_velocity_m_s
                        .iter()
                        .fold(0.0_f64, |s, u| s.max(u.abs())),
                );
            } else {
                close(field.pressure_pa, C64::ZERO, reference.pressure_pa.abs());
            }
        }
    }
}

#[test]
fn g1_rotating_wedge_images_enforce_rigid_and_release_faces() {
    for order in [3, 4] {
        let angle = core::f64::consts::PI / order as f64;
        let (sine, cosine) = angle.sin_cos();
        let n2 = [sine, -cosine, 0.0];
        for boundary in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let planes = [
                IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundary).unwrap(),
                // The point on face two is also our receiver: its distance is
                // exactly zero without snapping the admitted boundary.
                IdealAcousticPlane::new([cosine, sine, 0.0], n2, boundary).unwrap(),
            ];
            for receiver in [[0.7, 0.0, 0.2], [cosine, sine, 0.0]] {
                let source_position = [0.3 * (angle / 2.0).cos(), 0.3 * (angle / 2.0).sin(), -0.1];
                let direct = owner(51.0, receiver, &[source_position]);
                let bounded = direct.clone().with_wedge(planes, order, 2 * order).unwrap();
                let field = bounded.evaluate(&[source()]).unwrap();
                let reference = direct.evaluate(&[source()]).unwrap();
                if boundary == IdealPlaneBoundary::Rigid {
                    let normal = if receiver[1] == 0.0 {
                        [0.0, 1.0, 0.0]
                    } else {
                        n2
                    };
                    let velocity = (0..3).fold(C64::ZERO, |v, i| {
                        v + field.particle_velocity_m_s[i].scale(normal[i])
                    });
                    close(
                        velocity,
                        C64::ZERO,
                        reference
                            .particle_velocity_m_s
                            .iter()
                            .fold(0.0_f64, |s, u| s.max(u.abs())),
                    );
                } else {
                    close(field.pressure_pa, C64::ZERO, reference.pressure_pa.abs());
                }
                assert_eq!(bounded.reflection_wedge_order(), Some(order));
            }
        }
    }
}

#[test]
fn g4_reflection_geometry_budgets_and_cancellation_refuse_without_partial_results() {
    let plane =
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    let direct = owner(127.0, [0.4, 0.5, 0.6], &[[0.2, 0.3, 0.4]]);
    assert!(direct.clone().with_orthogonal_planes(&[plane], 1).is_err());
    assert!(direct.clone().with_orthogonal_planes(&[], 2).is_err());
    assert!(
        direct
            .clone()
            .with_orthogonal_planes(&[plane, plane], 4)
            .is_err()
    );
    assert!(
        owner(127.0, [0.4, -0.5, 0.6], &[[0.2, 0.3, 0.4]])
            .with_orthogonal_planes(&[plane], 2)
            .is_err()
    );
    assert!(
        owner(127.0, [0.4, 0.5, 0.6], &[[0.2, 0.0, 0.4]])
            .with_orthogonal_planes(&[plane], 2)
            .is_err()
    );
    assert!(
        direct
            .clone()
            .with_orthogonal_planes(&[plane], 2)
            .unwrap()
            .with_orthogonal_planes(&[plane], 2)
            .is_err()
    );
    for stop in [1, 5, 6] {
        let mut polls = 0;
        assert!(matches!(
            direct
                .clone()
                .with_orthogonal_planes_cancellable(&[plane], 2, || {
                    polls += 1;
                    polls == stop
                }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let again = direct.with_orthogonal_planes(&[plane], 2).unwrap();
    assert!(again.evaluate(&[source()]).is_ok());
}
