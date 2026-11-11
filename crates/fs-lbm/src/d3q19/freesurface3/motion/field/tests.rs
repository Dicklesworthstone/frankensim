use super::*;

fn tank(wall: PlanarWall3) -> FreeSurface3 {
    FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, y, z| {
            let coordinate = [x, y, z][wall.axis()];
            if wall.upper() {
                coordinate >= 4
            } else {
                coordinate < 4
            }
        },
    )
    .unwrap()
}

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 512,
        maximum_wall_faces: 64,
        maximum_relative_mass_defect: 1e-12,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn field() -> Vec<[f64; 2]> {
    (0..64)
        .map(|i| {
            let a = (i % 8) as f64 - 3.5;
            let b = (i / 8) as f64 - 3.5;
            [
                0.001 + 0.0004 * (a + 0.3 * b),
                -0.0015 + 0.0003 * (b - 0.2 * a),
            ]
        })
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn moments(surface: &FreeSurface3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if let Some(sample) = surface.sample_at(x, y, z) {
                    let mass = sample.fill_fraction * sample.density;
                    let momentum = sample.velocity.map(|u| mass * u);
                    let arm = [x, y, z].map(|p| p as f64 + 0.5);
                    let moment = cross3(core::array::from_fn(|d| arm[d] - origin[d]), momentum);
                    for d in 0..3 {
                        linear[d] += momentum[d];
                        angular[d] += moment[d];
                    }
                }
            }
        }
    }
    (linear, angular)
}

#[test]
fn g1_nonuniform_six_plane_shear_stretch_matches_independent_link_and_liquid_moments() {
    let velocities = field();
    let origin = [-1.7, 0.2, 4.3];
    let gate = CancelGate::new_clock_free();
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        let original = tank(wall);
        let motion = PlanarWallVelocityField3 {
            wall,
            moment_origin: origin,
            tangential_velocity: &velocities,
        };
        let (candidate, report) = original
            .planar_wall_velocity_field_candidate(motion, controls(), &gate)
            .unwrap();
        let mut stationary = original.clone();
        stationary.step().unwrap();
        let before = moments(&stationary, origin);
        let after = moments(&candidate, origin);
        let axis = wall.axis();
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let inward = if wall.upper() { -1.0 } else { 1.0 };
        let mut total_impulse = [0.0; 3];
        let mut total_angular = [0.0; 3];
        let mut total_work = 0.0;
        for (index, face) in report.faces.iter().enumerate() {
            let mut u = [0.0; 3];
            u[a] = velocities[index][0];
            u[b] = velocities[index][1];
            let mut expected_impulse = [0.0; 3];
            let mut expected_angular = [0.0; 3];
            let mut expected_couple = [0.0; 3];
            // Independent D3Q19 incoming half-space: one axial weight 1/18,
            // and the four axial/tangent diagonal weights 1/36.
            for (tangent, sign) in [(axis, 0.0), (a, -1.0), (a, 1.0), (b, -1.0), (b, 1.0)] {
                let mut e = [0.0; 3];
                e[axis] = inward;
                if tangent != axis {
                    e[tangent] = sign;
                }
                let w = if tangent == axis {
                    1.0 / 18.0
                } else {
                    1.0 / 36.0
                };
                let strength = 2.0 * w + 6.0 * w * (e[a] * u[a] + e[b] * u[b]);
                let impulse = e.map(|value| value * strength);
                let point: [f64; 3] =
                    core::array::from_fn(|d| face.liquid_cell[d] as f64 + 0.5 - 0.5 * e[d]);
                let angular = cross3(core::array::from_fn(|d| point[d] - origin[d]), impulse);
                let couple = cross3(core::array::from_fn(|d| point[d] - face.center[d]), impulse);
                for d in 0..3 {
                    expected_impulse[d] += impulse[d];
                    expected_angular[d] += angular[d];
                    expected_couple[d] += couple[d];
                }
            }
            for d in 0..3 {
                close(face.fluid_impulse[d], expected_impulse[d]);
                close(face.fluid_angular_impulse[d], expected_angular[d]);
                close(face.fluid_couple_at_face[d], expected_couple[d]);
                assert_eq!(face.wall_impulse[d], -face.fluid_impulse[d]);
                assert_eq!(face.wall_angular_impulse[d], -face.fluid_angular_impulse[d]);
                assert_eq!(face.wall_couple_at_face[d], -face.fluid_couple_at_face[d]);
                total_impulse[d] += expected_impulse[d];
                total_angular[d] += expected_angular[d];
            }
            let force_moment = cross3(
                core::array::from_fn(|d| face.center[d] - origin[d]),
                face.fluid_impulse,
            );
            for d in 0..3 {
                close(
                    force_moment[d] + face.fluid_couple_at_face[d],
                    face.fluid_angular_impulse[d],
                );
            }
            close(face.fluid_tangential_response[0], 1.0 / 3.0);
            close(face.fluid_tangential_response[1], 1.0 / 3.0);
            let work = (u[a] * u[a] + u[b] * u[b]) / 3.0;
            close(face.fluid_work, work);
            assert_eq!(face.wall_work, -face.fluid_work);
            total_work += work;
            assert!(face.relative_mass_defect <= 1e-12);
            let sample = candidate
                .sample_at(
                    face.liquid_cell[0],
                    face.liquid_cell[1],
                    face.liquid_cell[2],
                )
                .unwrap();
            close(sample.velocity[a], u[a] / 3.0);
            close(sample.velocity[b], u[b] / 3.0);
        }
        for d in 0..3 {
            close(report.fluid_impulse[d], total_impulse[d]);
            close(report.fluid_angular_impulse[d], total_angular[d]);
            close(
                after.0[d] - before.0[d],
                report.fluid_impulse[d] - report.reference_pressure_fluid_impulse[d],
            );
            close(
                after.1[d] - before.1[d],
                report.fluid_angular_impulse[d]
                    - report.reference_pressure_fluid_angular_impulse[d],
            );
        }
        close(report.fluid_work, total_work);
        close(report.mass_after, report.mass_before);
        assert!(
            report
                .faces
                .iter()
                .any(|f| f.fluid_couple_at_face.iter().any(|v| v.abs() > 1e-6))
        );
        assert_eq!(original.grid.state_bits(), tank(wall).grid.state_bits());
    }
}

#[test]
fn g3_face_response_replays_current_populations_and_moments_shift_with_the_origin() {
    let wall = PlanarWall3::ZMin;
    let gate = CancelGate::new_clock_free();
    let mut surface = tank(wall);
    let velocities = field();
    let motion = PlanarWallVelocityField3 {
        wall,
        moment_origin: [0.0; 3],
        tangential_velocity: &velocities,
    };
    surface
        .step_with_planar_wall_velocity_field(motion, controls(), &gate)
        .unwrap();
    let (_, pilot) = surface
        .planar_wall_velocity_field_candidate(motion, controls(), &gate)
        .unwrap();
    let mut changed = velocities.clone();
    changed[13][0] += 0.0007;
    changed[13][1] -= 0.0003;
    let second_motion = PlanarWallVelocityField3 {
        tangential_velocity: &changed,
        ..motion
    };
    let (after, second) = surface
        .planar_wall_velocity_field_candidate(second_motion, controls(), &gate)
        .unwrap();
    for (index, (first, final_face)) in pilot.faces.iter().zip(&second.faces).enumerate() {
        let rho = surface
            .sample_at(
                first.liquid_cell[0],
                first.liquid_cell[1],
                first.liquid_cell[2],
            )
            .unwrap()
            .density;
        for d in 0..2 {
            close(first.fluid_tangential_response[d], rho / 3.0);
            close(
                final_face.fluid_impulse[d],
                first.fluid_impulse[d]
                    + first.fluid_tangential_response[d]
                        * (changed[index][d] - velocities[index][d]),
            );
        }
        close(final_face.fluid_impulse[2], first.fluid_impulse[2]);
        let actual_work = final_face
            .velocity
            .iter()
            .zip(final_face.fluid_impulse)
            .map(|(u, j)| u * j)
            .sum();
        close(final_face.fluid_work, actual_work);
    }
    let shift = [2.3, -1.1, 0.7];
    let (shifted, moved) = surface
        .planar_wall_velocity_field_candidate(
            PlanarWallVelocityField3 {
                moment_origin: shift,
                ..second_motion
            },
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(after.grid.state_bits(), shifted.grid.state_bits());
    let correction = cross3(shift, second.fluid_impulse);
    for d in 0..3 {
        close(
            moved.fluid_angular_impulse[d],
            second.fluid_angular_impulse[d] - correction[d],
        );
    }
    for (a, b) in second.faces.iter().zip(&moved.faces) {
        assert_eq!(a.fluid_couple_at_face, b.fluid_couple_at_face);
    }
    let mut continuation = after.clone();
    let mut replay = after;
    for _ in 0..3 {
        let a = continuation
            .step_with_planar_wall_velocity_field(second_motion, controls(), &gate)
            .unwrap();
        let b = replay
            .step_with_planar_wall_velocity_field(second_motion, controls(), &gate)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(continuation.grid.state_bits(), replay.grid.state_bits());
    }
}

#[test]
fn g3_uniform_and_zero_fields_preserve_original_native_population_arithmetic() {
    let gate = CancelGate::new_clock_free();
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        for tangent in [[0.0; 2], [0.003, -0.002]] {
            let velocities = vec![tangent; 64];
            let mut field = tank(wall);
            let mut original = field.clone();
            let mut velocity = [0.0; 3];
            velocity[(wall.axis() + 1) % 3] = tangent[0];
            velocity[(wall.axis() + 2) % 3] = tangent[1];
            let old = original
                .step_with_planar_wall_motion(wall, velocity, 1e-12, &gate)
                .unwrap();
            let new = field
                .step_with_planar_wall_velocity_field(
                    PlanarWallVelocityField3 {
                        wall,
                        moment_origin: [0.0; 3],
                        tangential_velocity: &velocities,
                    },
                    controls(),
                    &gate,
                )
                .unwrap();
            assert_eq!(field.grid.state_bits(), original.grid.state_bits());
            assert_eq!(field.mass_snapshot(), original.mass_snapshot());
            assert_eq!(field.conversions(), original.conversions());
            for d in 0..3 {
                close(new.fluid_impulse[d], old.fluid_impulse[d]);
            }
            close(new.fluid_work, old.fluid_work);
        }
    }
}

#[test]
fn g4_field_shape_resource_cancel_and_late_moment_refusals_preserve_exact_retry() {
    let gate = CancelGate::new_clock_free();
    let mut surface = tank(PlanarWall3::ZMin);
    let velocities = field();
    let motion = PlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin: [0.0; 3],
        tangential_velocity: &velocities,
    };
    let before = surface.grid.state_bits();
    let masses = surface.mass_snapshot();
    let counters = surface.conversions();
    for bad in [
        PlanarWallFieldControls3 {
            maximum_domain_cells: 511,
            ..controls()
        },
        PlanarWallFieldControls3 {
            maximum_wall_faces: 63,
            ..controls()
        },
        PlanarWallFieldControls3 {
            maximum_relative_face_mass_defect: f64::NAN,
            ..controls()
        },
    ] {
        assert!(
            surface
                .step_with_planar_wall_velocity_field(motion, bad, &gate)
                .is_err()
        );
    }
    assert!(
        surface
            .step_with_planar_wall_velocity_field(
                PlanarWallVelocityField3 {
                    tangential_velocity: &velocities[..63],
                    ..motion
                },
                controls(),
                &gate,
            )
            .is_err()
    );
    let mut fast = velocities.clone();
    // The wall field is refused at `speed / sqrt(1/3) >= MACH_LIMIT` (field.rs, 9c14315625, 2026-11-03), i.e. from
    // 0.3 x sqrt(1/3) = 0.173; the old literal 0.1 is below it. Twice the limit follows it if it ever changes.
    fast[31] = [2.0 * MACH_LIMIT * CS2_3.sqrt(), 0.0];
    assert!(
        surface
            .step_with_planar_wall_velocity_field(
                PlanarWallVelocityField3 {
                    tangential_velocity: &fast,
                    ..motion
                },
                controls(),
                &gate,
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        surface.step_with_planar_wall_velocity_field(motion, controls(), &cancelled),
        Err(FreeSurfaceError3::Cancelled)
    ));
    // A finite origin passes field admission; actual diagonal-link cross
    // products overflow only after the original candidate has begun streaming.
    assert!(
        surface
            .step_with_planar_wall_velocity_field(
                PlanarWallVelocityField3 {
                    moment_origin: [f64::MAX; 3],
                    ..motion
                },
                controls(),
                &gate,
            )
            .is_err()
    );
    assert_eq!(surface.grid.state_bits(), before);
    assert_eq!(surface.mass_snapshot(), masses);
    assert_eq!(surface.conversions(), counters);
    let mut direct = tank(PlanarWall3::ZMin);
    let a = surface
        .step_with_planar_wall_velocity_field(motion, controls(), &gate)
        .unwrap();
    let b = direct
        .step_with_planar_wall_velocity_field(motion, controls(), &gate)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(surface.grid.state_bits(), direct.grid.state_bits());
}
