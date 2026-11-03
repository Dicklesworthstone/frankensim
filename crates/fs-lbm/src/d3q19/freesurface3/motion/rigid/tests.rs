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

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-12 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn base() -> PlanarRigidWallMotion3 {
    PlanarRigidWallMotion3 {
        wall: PlanarWall3::ZMin,
        origin: [4.0, 4.0, 0.0],
        velocity_at_origin: [0.0; 3],
        angular_speed: 0.002,
    }
}

fn angular_momentum(surface: &FreeSurface3, origin: [f64; 3]) -> [f64; 3] {
    let mut result = [0.0; 3];
    for (&key, tile) in &surface.tiles {
        let slot = surface.grid.slot_of(key).unwrap();
        for (lane, &cell) in tile.cells.iter().enumerate() {
            if cell == Cell3::Gas {
                continue;
            }
            let (x, y, z) = FreeSurface3::global_of(key, lane);
            let sample = surface
                .sample_at(x as usize, y as usize, z as usize)
                .unwrap();
            let mass: f64 = if cell == Cell3::Fluid {
                surface.grid.populations(slot, lane).iter().sum()
            } else {
                tile.mass[lane]
            };
            let position = [x, y, z].map(|coordinate| coordinate as f64 + 0.5);
            let moment = cross(
                core::array::from_fn(|d| position[d] - origin[d]),
                sample.velocity.map(|v| v * mass),
            );
            for d in 0..3 {
                result[d] += moment[d];
            }
        }
    }
    result
}

#[test]
fn uniform_rigid_entry_preserves_all_six_original_wall_population_paths() {
    let gate = CancelGate::new_clock_free();
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        let initial = tank(wall);
        for speed in [0.0, 0.01, -0.01] {
            let mut velocity = [0.0; 3];
            velocity[(wall.axis() + 1) % 3] = speed;
            let mut uniform = initial.clone();
            let old = uniform
                .step_with_planar_wall_motion(wall, velocity, 1e-12, &gate)
                .unwrap();
            let (rigid, new) = initial
                .planar_rigid_wall_motion_candidate(
                    PlanarRigidWallMotion3 {
                        wall,
                        origin: [2.0, 3.0, 4.0],
                        velocity_at_origin: velocity,
                        angular_speed: 0.0,
                    },
                    1e-12,
                    &gate,
                )
                .unwrap();
            assert_eq!(rigid.grid.state_bits(), uniform.grid.state_bits());
            assert_eq!(rigid.mass_snapshot(), uniform.mass_snapshot());
            assert_eq!(rigid.conversions(), uniform.conversions());
            assert_eq!(new.fluid_impulse, old.fluid_impulse);
            assert_eq!(new.wall_impulse, old.wall_impulse);
            assert_eq!(
                new.reference_pressure_fluid_impulse,
                old.reference_pressure_fluid_impulse
            );
            near(new.fluid_work, old.fluid_work);
            assert_eq!(new.mass_after, old.mass_after);
            assert_eq!(new.maximum_speed, old.maximum_speed);
        }
    }
}

#[test]
fn rotation_matches_independent_d3q19_impedance_and_population_angular_momentum() {
    let gate = CancelGate::new_clock_free();
    let initial = tank(PlanarWall3::ZMin);
    let motion = base();
    let before = angular_momentum(&initial, motion.origin);
    let (after, report) = initial
        .planar_rigid_wall_motion_candidate(motion, 1e-12, &gate)
        .unwrap();
    // Sum receiving-cell radii, independently of the wall's response code.
    // Each tangent has the original D3Q19 impedance rho/3 per full wall cell.
    let mut impedance = 0.0;
    for y in 0..8 {
        for x in 0..8 {
            impedance += ((x as f64 + 0.5 - 4.0).powi(2) + (y as f64 + 0.5 - 4.0).powi(2)) / 3.0;
        }
    }
    near(impedance, 224.0);
    near(report.fluid_response[2][2], impedance);
    near(
        report.fluid_angular_impulse[2],
        impedance * motion.angular_speed,
    );
    near(
        angular_momentum(&after, motion.origin)[2] - before[2],
        report.fluid_angular_impulse[2],
    );
    near(
        report.fluid_work,
        motion.angular_speed * report.fluid_angular_impulse[2],
    );
    near(report.mass_before, report.mass_after);
    for i in 0..3 {
        assert_eq!(
            report.wall_angular_impulse[i],
            -report.fluid_angular_impulse[i]
        );
        for j in 0..3 {
            assert_eq!(report.fluid_response[i][j], report.fluid_response[j][i]);
        }
    }
    near(report.fluid_response[0][0], 64.0 / 3.0);
    near(report.fluid_response[1][1], 64.0 / 3.0);
    near(report.fluid_response[0][2], 0.0);
    near(report.fluid_response[1][2], 0.0);
    let changed = PlanarRigidWallMotion3 {
        velocity_at_origin: [0.005, -0.003, 0.0],
        angular_speed: -0.001,
        ..motion
    };
    let (_, replay) = initial
        .planar_rigid_wall_motion_candidate(changed, 1e-12, &gate)
        .unwrap();
    let delta = [
        changed.velocity_at_origin[0],
        changed.velocity_at_origin[1],
        changed.angular_speed - motion.angular_speed,
    ];
    let old = [
        report.fluid_impulse[0],
        report.fluid_impulse[1],
        report.fluid_angular_impulse[2],
    ];
    let actual = [
        replay.fluid_impulse[0],
        replay.fluid_impulse[1],
        replay.fluid_angular_impulse[2],
    ];
    for i in 0..3 {
        near(
            actual[i],
            old[i]
                + (0..3)
                    .map(|j| report.fluid_response[i][j] * delta[j])
                    .sum::<f64>(),
        );
    }
}

#[test]
fn origin_shift_preserves_local_motion_work_and_transforms_both_raw_and_reference_moments() {
    let gate = CancelGate::new_clock_free();
    let initial = tank(PlanarWall3::ZMin);
    let motion = PlanarRigidWallMotion3 {
        velocity_at_origin: [1.0 / 256.0, -1.0 / 512.0, 0.0],
        angular_speed: 1.0 / 1024.0,
        ..base()
    };
    let shift = [1.0, -2.0, 3.0];
    let shifted = PlanarRigidWallMotion3 {
        origin: core::array::from_fn(|d| motion.origin[d] + shift[d]),
        velocity_at_origin: [
            motion.velocity_at_origin[0] - motion.angular_speed * shift[1],
            motion.velocity_at_origin[1] + motion.angular_speed * shift[0],
            0.0,
        ],
        ..motion
    };
    let (a, old) = initial
        .planar_rigid_wall_motion_candidate(motion, 1e-12, &gate)
        .unwrap();
    let (b, new) = initial
        .planar_rigid_wall_motion_candidate(shifted, 1e-12, &gate)
        .unwrap();
    assert_eq!(a.grid.state_bits(), b.grid.state_bits());
    assert_eq!(old.fluid_impulse, new.fluid_impulse);
    near(new.fluid_work, old.fluid_work);
    let shifted_moment = cross(shift, old.fluid_impulse);
    let shifted_reference = cross(shift, old.reference_pressure_fluid_impulse);
    for d in 0..3 {
        near(
            new.fluid_angular_impulse[d],
            old.fluid_angular_impulse[d] - shifted_moment[d],
        );
        near(
            new.reference_pressure_fluid_angular_impulse[d],
            old.reference_pressure_fluid_angular_impulse[d] - shifted_reference[d],
        );
    }
    near(
        new.fluid_work,
        (0..3)
            .map(|d| shifted.velocity_at_origin[d] * new.fluid_impulse[d])
            .sum::<f64>()
            + shifted.angular_speed * new.fluid_angular_impulse[2],
    );
}

#[test]
fn local_speed_geometry_and_cancellation_refuse_without_publication() {
    let mut initial = tank(PlanarWall3::ZMin);
    let before = initial.grid.state_bits();
    let gate = CancelGate::new_clock_free();
    for motion in [
        PlanarRigidWallMotion3 {
            angular_speed: 1.0,
            ..base()
        },
        PlanarRigidWallMotion3 {
            velocity_at_origin: [0.0, 0.0, 0.01],
            ..base()
        },
        PlanarRigidWallMotion3 {
            origin: [f64::INFINITY, 0.0, 0.0],
            ..base()
        },
        PlanarRigidWallMotion3 {
            wall: PlanarWall3::ZMax,
            ..base()
        },
    ] {
        assert!(
            initial
                .step_with_planar_rigid_wall_motion(motion, 1e-12, &gate)
                .is_err()
        );
        assert_eq!(initial.grid.state_bits(), before);
    }
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        initial.step_with_planar_rigid_wall_motion(base(), 1e-12, &stopped),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert_eq!(initial.grid.state_bits(), before);
}
