//! Prescribed rigid-frame motion through original D3Q19 collision and streaming.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, CollisionModel3, E3, FaceBoundary3, Q3,
    RigidFrame3, collide_cell3, equilibrium3,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn closed() -> BoundarySpec3 {
    BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6])
}
fn moments(f: [f64; Q3], force: [f64; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let directions = E3.map(|(a, b, c)| [f64::from(a), f64::from(b), f64::from(c)]);
    let rho: f64 = f.iter().sum();
    let momentum: [f64; 3] = core::array::from_fn(|a| {
        f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]
    });
    let velocity = momentum.map(|p| p / rho);
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
        })
    });
    (rho, momentum, stress)
}
fn local(grid: &BoundaryGrid3, [x, y, z]: [usize; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    moments(grid.populations(x, y, z), grid.force_density(x, y, z))
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<([f64; Q3], [f64; 3])> {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                result.push((grid.populations(x, y, z), grid.force_density(x, y, z)));
            }
        }
    }
    result
}
fn total_momentum(grid: &BoundaryGrid3, origin: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (mut linear, mut angular) = ([0.0; 3], [0.0; 3]);
    let [nx, ny, nz] = grid.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let momentum = local(grid, [x, y, z]).1;
                let arm = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let moment = cross(arm, momentum);
                for a in 0..3 {
                    linear[a] += momentum[a];
                    angular[a] += moment[a];
                }
            }
        }
    }
    (linear, angular)
}

#[test]
fn arbitrary_axis_coriolis_rotates_without_work_in_all_original_collision_models() {
    let models = [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.1,
            higher_order_rate: 1.4,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.1,
            third_order_rate: 1.3,
            fourth_order_rate: 1.5,
        },
    ];
    let gate = CancelGate::new();
    let omega = [0.04, -0.03, 0.02];
    let frame = RigidFrame3 {
        origin: [1.5, 2.5, 1.5],
        angular_velocity: omega,
        ..Default::default()
    };
    for model in models {
        let mut grid = BoundaryGrid3::with_collision_model(4, 4, 4, model, [0.0; 3], closed());
        grid.initialize_uniform_equilibrium(1.3, [0.03, -0.015, 0.021], &gate)
            .unwrap();
        let before = local(&grid, [1, 2, 1]);
        let velocity = before.1.map(|p| p / before.0);
        frame.prepare_step(&mut grid, [0.0; 3], &gate).unwrap();
        let force = grid.force_density(1, 2, 1);
        let post = collide_cell3(grid.populations(1, 2, 1), model, force).unwrap();
        let after = moments(post, force);
        let next_velocity = after.1.map(|p| p / after.0);
        // Independent Rodrigues oracle: midpoint Coriolis turns by -2 atan(|omega|).
        let magnitude = dot(omega, omega).sqrt();
        let axis = omega.map(|w| w / magnitude);
        let angle = 2.0 * magnitude.atan();
        let parallel = dot(axis, velocity);
        let turn = cross(axis, velocity);
        for a in 0..3 {
            let expected = axis[a] * parallel + angle.cos() * (velocity[a] - axis[a] * parallel)
                - angle.sin() * turn[a];
            near(next_velocity[a], expected, 8e-15);
        }
        near(
            dot(next_velocity, next_velocity),
            dot(velocity, velocity),
            2e-15,
        );
        near(dot(axis, next_velocity), parallel, 8e-15);
        let midpoint = core::array::from_fn(|a| 0.5 * (velocity[a] + next_velocity[a]));
        near(dot(force, midpoint), 0.0, 2e-15);
        assert!(dot(force, force) > 1e-7);
    }
}

#[test]
fn euler_centrifugal_and_origin_acceleration_use_native_half_cell_positions() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], closed());
    let frame = RigidFrame3 {
        origin: [1.5, 1.5, 2.5],
        angular_velocity: [0.01, -0.02, 0.03],
        angular_acceleration: [0.002, 0.003, -0.001],
        origin_acceleration: [0.001, -0.002, 0.005],
    };
    // At centre (2.5,3.5,1.5), r=(1,2,-1). Centrifugal=(.002,.0016,.0004),
    // Euler=(.001,-.001,-.001), and g-a_origin=(.003,-.004,-.002).
    let expected = [0.006, -0.0034, -0.0026];
    let velocity = expected.map(|a| -0.5 * a); // u_mid=0, so no Coriolis contribution.
    grid.set_populations(2, 3, 1, &equilibrium3(1.25, velocity));
    frame
        .prepare_step(&mut grid, [0.004, -0.006, 0.003], &gate)
        .unwrap();
    for a in 0..3 {
        near(grid.force_density(2, 3, 1)[a], 1.25 * expected[a], 2e-15);
        near(grid.velocity(2, 3, 1)[a], velocity[a], 2e-15);
    }
}

#[test]
fn closed_box_streaming_and_wall_tractions_close_linear_and_angular_impulse() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(8, 8, 8, 0.9, [0.0; 3], closed());
    grid.initialize_uniform_equilibrium(1.2, [0.004, -0.003, 0.002], &gate)
        .unwrap();
    let frame = RigidFrame3 {
        origin: [4.0; 3],
        angular_velocity: [0.003, -0.002, 0.004],
        angular_acceleration: [2e-5, -3e-5, 4e-5],
        origin_acceleration: [2e-5, 1e-5, -1e-5],
    };
    frame
        .prepare_step(&mut grid, [1e-5, -4e-5, 2e-5], &gate)
        .unwrap();
    let origin = [2.7, -1.2, 0.4];
    let (before, angular_before) = total_momentum(&grid, origin);
    let initial_mass = grid.total_mass();
    let (mut applied, mut applied_angular) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let force = grid.force_density(x, y, z);
                let arm = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let torque = cross(arm, force);
                for a in 0..3 {
                    applied[a] += force[a];
                    applied_angular[a] += torque[a];
                }
            }
        }
    }
    let walls = grid.boundary_links();
    let exchange = grid.step_with_wall_momentum(&walls, origin);
    let (after, angular_after) = total_momentum(&grid, origin);
    for a in 0..3 {
        near(
            after[a] - before[a] + exchange.wall_impulse[a],
            applied[a],
            2e-12,
        );
        near(
            angular_after[a] - angular_before[a] + exchange.wall_angular_impulse[a],
            applied_angular[a],
            2e-11,
        );
    }
    near(grid.total_mass(), initial_mass, 2e-10);
    assert!(dot(applied_angular, applied_angular) > 0.01);
    assert!(exchange.measured_links > 0);
}

#[test]
fn force_replacement_preserves_physical_momentum_and_nonequilibrium_stress() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [1e-4, -2e-4, 3e-4], closed());
    grid.set_force_field(|[x, y, z]| [x as f64 * 1e-5, y as f64 * 2e-5, -(z as f64) * 1e-5]);
    grid.set_force_field_scale(1.5);
    grid.initialize_uniform_equilibrium(1.3, [0.03, -0.02, 0.01], &gate)
        .unwrap();
    let mut f = grid.populations(2, 1, 3);
    f[1] += 0.002;
    f[2] += 0.002;
    f[0] -= 0.004;
    grid.set_populations(2, 1, 3, &f);
    let original = local(&grid, [2, 1, 3]);
    let frame = RigidFrame3 {
        origin: [2.0; 3],
        angular_velocity: [0.003, -0.002, 0.004],
        angular_acceleration: [1e-5, 2e-5, -1e-5],
        origin_acceleration: [2e-5, 0.0, -1e-5],
    };
    frame
        .prepare_step(&mut grid, [0.0, -1e-4, 0.0], &gate)
        .unwrap();
    assert_eq!(grid.uniform_force_density(), [0.0; 3]);
    let selected_force = grid.force_density(2, 1, 3);
    // An older selected shape must not survive the frame owner's replacement.
    grid.set_force_field_scale(2.0);
    for a in 0..3 {
        near(grid.force_density(2, 1, 3)[a], 2.0 * selected_force[a], 0.0);
    }
    for motion in [frame, RigidFrame3::default()] {
        motion.prepare_step(&mut grid, [0.0; 3], &gate).unwrap();
        let next = local(&grid, [2, 1, 3]);
        near(next.0, original.0, 4e-15);
        for a in 0..3 {
            near(next.1[a], original.1[a], 4e-15);
            for b in 0..3 {
                near(next.2[a][b], original.2[a][b], 4e-15);
            }
        }
    }
    assert_eq!(grid.force_density(2, 1, 3), [0.0; 3]);
}

#[test]
fn periodicity_requires_translation_invariant_inertial_forces() {
    let gate = CancelGate::new();
    let spin = RigidFrame3 {
        origin: [2.0; 3],
        angular_velocity: [0.0, 0.0, 0.002],
        angular_acceleration: [0.0, 0.0, 1e-5],
        ..Default::default()
    };
    let mut torus = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    let before = snapshot(&torus);
    assert!(matches!(
        spin.prepare_step(&mut torus, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Invalid(_))
    ));
    assert_eq!(snapshot(&torus), before);
    RigidFrame3::default()
        .prepare_step(&mut torus, [1e-5, -2e-5, 3e-5], &gate)
        .unwrap();
    let wall = FaceBoundary3::stationary_wall();
    let mut axial = BoundaryGrid3::new(
        4,
        4,
        4,
        0.9,
        [0.0; 3],
        BoundarySpec3::new([
            wall,
            wall,
            wall,
            wall,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
        ]),
    );
    spin.prepare_step(&mut axial, [0.0; 3], &gate).unwrap();
    assert_eq!(axial.force_density(1, 2, 0), axial.force_density(1, 2, 3));
    let incompatible_alpha = RigidFrame3 {
        angular_acceleration: [1e-5, 0.0, 0.0],
        ..spin
    };
    assert!(
        incompatible_alpha
            .prepare_step(&mut axial, [0.0; 3], &gate)
            .is_err()
    );
    let mut open = BoundaryGrid3::new(
        4,
        4,
        4,
        0.9,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(1.0, 1.0),
    );
    assert!(
        RigidFrame3::default()
            .prepare_step(&mut open, [0.0; 3], &gate)
            .is_err()
    );
}

#[test]
fn late_unphysical_cell_invalid_motion_and_cancellation_publish_nothing() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(8, 8, 8, 0.9, [1e-4, -2e-4, 3e-4], closed());
    grid.set_force_field(|[x, y, z]| [x as f64 * 1e-5, y as f64 * 2e-5, -(z as f64) * 1e-5]);
    grid.set_force_field_scale(1.5);
    grid.set_populations(7, 7, 7, &[-1.0; Q3]);
    let before = snapshot(&grid);
    let frame = RigidFrame3 {
        origin: [4.0; 3],
        angular_velocity: [0.003, -0.002, 0.001],
        ..Default::default()
    };
    assert_eq!(
        frame.prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Unphysical { cell: [7, 7, 7] })
    );
    assert_eq!(snapshot(&grid), before);
    assert_eq!(grid.uniform_force_density(), [1e-4, -2e-4, 3e-4]);
    let overflow = RigidFrame3 {
        angular_velocity: [f64::MAX, 0.0, 0.0],
        ..frame
    };
    assert!(matches!(
        overflow.prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Invalid(_))
    ));
    gate.request();
    assert_eq!(
        frame.prepare_step(&mut grid, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Cancelled)
    );
    assert_eq!(snapshot(&grid), before);
}
