//! Native finite-inertia disk response, momentum balance, and atomic refusal.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, BoundaryStepError3, CollisionModel3, E3, Face3, FaceBoundary3,
    Q3, WallRotor3, WallRotorError3, equilibrium3,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (error {:e})",
        actual - expected
    );
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
fn axis_and_min(face: Face3) -> (usize, bool) {
    match face {
        Face3::XMin => (0, true),
        Face3::XMax => (0, false),
        Face3::YMin => (1, true),
        Face3::YMax => (1, false),
        Face3::ZMin => (2, true),
        Face3::ZMax => (2, false),
    }
}
fn center(n: usize, face: Face3) -> [f64; 3] {
    let (axis, min) = axis_and_min(face);
    let mut center = [n as f64 / 2.0; 3];
    center[axis] = if min { 0.0 } else { n as f64 };
    center
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !grid.is_solid(x, y, z) {
                    result.push(grid.populations(x, y, z).map(f64::to_bits));
                }
            }
        }
    }
    result
}
fn inventory(grid: &BoundaryGrid3, origin: [f64; 3]) -> (f64, [f64; 3], [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    let mut mass = 0.0;
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let f = grid.populations(x, y, z);
                let force = grid.force_density(x, y, z);
                let mut momentum = force.map(|f| 0.5 * f);
                for q in 0..Q3 {
                    let e = [E3[q].0, E3[q].1, E3[q].2];
                    mass += f[q];
                    for i in 0..3 {
                        momentum[i] += f64::from(e[i]) * f[q];
                    }
                }
                let cell = [x, y, z];
                let moment = cross(
                    core::array::from_fn(|i| cell[i] as f64 + 0.5 - origin[i]),
                    momentum,
                );
                for i in 0..3 {
                    linear[i] += momentum[i];
                    angular[i] += moment[i];
                }
            }
        }
    }
    (mass, linear, angular)
}
// Geometric second area moment, independently enumerated in face coordinates.
// The density is one in the resting-fluid experiment.
fn resting_damping(n: usize, face: Face3, origin: [f64; 3], radius: f64) -> (usize, f64) {
    let (axis, _) = axis_and_min(face);
    let tangent: Vec<_> = (0..3).filter(|&a| a != axis).collect();
    let mut count = 0;
    let mut second_moment = 0.0;
    for j in 0..n {
        for i in 0..n {
            let a = i as f64 + 0.5 - origin[tangent[0]];
            let b = j as f64 + 0.5 - origin[tangent[1]];
            if a * a + b * b <= radius * radius {
                count += 1;
                second_moment += a * a + b * b;
            }
        }
    }
    (count, second_moment / 3.0)
}

#[test]
fn native_link_slope_and_midpoint_spin_match_all_faces_and_two_radii() {
    let n = 8;
    let gate = CancelGate::new();
    for face in Face3::ALL {
        for radius in [1.4, 2.6] {
            let origin = center(n, face);
            let (count, damping) = resting_damping(n, face, origin, radius);
            let mut responses = Vec::new();
            for spin in [-0.003, 0.004] {
                let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
                let inertia = 37.0;
                let external = 0.015;
                let mut rotor =
                    WallRotor3::new(&grid, face, origin, radius, inertia, spin, n * n, &gate)
                        .unwrap();
                assert_eq!(rotor.patch_cell_count(), count);
                let receipt = grid
                    .step_with_finite_inertia_rotor(&mut rotor, external, origin, &gate)
                    .unwrap();
                let expected =
                    ((inertia - 0.5 * damping) * spin + external) / (inertia + 0.5 * damping);
                near(rotor.angular_velocity(), expected, 2e-16);
                near(
                    receipt.angular_velocity_midpoint,
                    0.5 * (spin + expected),
                    2e-16,
                );
                near(rotor.angle(), receipt.angular_velocity_midpoint, 0.0);
                near(
                    receipt.hydrodynamic_angular_impulse,
                    -damping * receipt.angular_velocity_midpoint,
                    2e-14,
                );
                near(
                    inertia * (expected - spin),
                    external + receipt.hydrodynamic_angular_impulse,
                    2e-14,
                );
                near(receipt.hydrodynamic.fluid_mass_change, 0.0, 1e-15);
                assert_eq!(receipt.hydrodynamic.measured_links, count * 5);
                near(receipt.mechanical_work_residual, 0.0, 2e-16);
                responses.push(receipt);
            }
            // Slope is measured from actual native moving-wall population exchange,
            // not from the scalar response's internal damping coefficient.
            let slope = (responses[1].hydrodynamic_angular_impulse
                - responses[0].hydrodynamic_angular_impulse)
                / (responses[1].angular_velocity_midpoint - responses[0].angular_velocity_midpoint);
            near(slope, -damping, 2e-11);
        }
    }
}

#[test]
fn small_rotor_inertia_uses_stable_midpoint_reaction() {
    let n = 8;
    let origin = center(n, Face3::ZMin);
    let (_, damping) = resting_damping(n, Face3::ZMin, origin, 2.6);
    let inertia = damping * 0.02;
    let spin = 0.004;
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
    let mut rotor =
        WallRotor3::new(&grid, Face3::ZMin, origin, 2.6, inertia, spin, n * n, &gate).unwrap();
    let receipt = grid
        .step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &gate)
        .unwrap();
    near(
        rotor.angular_velocity(),
        spin * (inertia - 0.5 * damping) / (inertia + 0.5 * damping),
        2e-16,
    );
    assert!(rotor.angular_velocity().abs() < spin);
    assert!(receipt.kinetic_energy_change < 0.0);
    near(receipt.mechanical_work_residual, 0.0, 1e-17);
}

#[test]
fn actual_fluid_rotor_housing_and_source_close_about_an_arbitrary_origin() {
    let gate = CancelGate::new();
    let models = [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::Trt {
            tau: 0.8,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 1.1,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 1.1,
            fourth_order_rate: 1.0,
        },
    ];
    let origin = [-1.7, 2.3, 9.4];
    let center = [4.25, 3.75, 0.0];
    let force = [2e-6, -3e-6, 1e-6];
    for model in models {
        let n = 8;
        let mut grid = BoundaryGrid3::with_collision_model(n, n, n, model, force, closed());
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let u = [
                        -0.002 * (y as f64 + 0.5 - center[1]) + 0.001,
                        0.002 * (x as f64 + 0.5 - center[0]),
                        0.0003 * (x as f64 - y as f64),
                    ];
                    let rho = 1.0 + 0.0001 * (x + y + z) as f64;
                    grid.set_populations(x, y, z, &equilibrium3(rho, u));
                }
            }
        }
        let inertia = 45.0;
        let mut rotor = WallRotor3::new(
            &grid,
            Face3::ZMin,
            center,
            2.1,
            inertia,
            -0.001,
            n * n,
            &gate,
        )
        .unwrap();
        for _ in 0..5 {
            let before = inventory(&grid, origin);
            let spin_before = rotor.angular_velocity();
            let external = 0.003;
            let receipt = grid
                .step_with_finite_inertia_rotor(&mut rotor, external, origin, &gate)
                .unwrap();
            let after = inventory(&grid, origin);
            near(after.0, before.0, 3e-10);
            let rotor_impulse = inertia * (rotor.angular_velocity() - spin_before);
            near(
                rotor_impulse,
                external + receipt.hydrodynamic_angular_impulse,
                1e-12,
            );
            for i in 0..3 {
                near(
                    after.1[i] - before.1[i] + receipt.stationary_walls.wall_impulse[i]
                        - receipt.bearing_impulse[i],
                    receipt.body_force_impulse[i],
                    2e-12,
                );
                near(
                    after.2[i] - before.2[i]
                        + if i == 2 { rotor_impulse } else { 0.0 }
                        + receipt.stationary_walls.wall_angular_impulse[i]
                        - receipt.bearing_angular_impulse[i],
                    receipt.body_force_angular_impulse[i] + if i == 2 { external } else { 0.0 },
                    6e-12,
                );
                near(
                    receipt.hydrodynamic.wall_velocity_mass_impulse[i],
                    0.0,
                    1e-16,
                );
                near(
                    receipt.hydrodynamic.wall_velocity_mass_angular_impulse[i],
                    0.0,
                    1e-15,
                );
            }
            let bearing_force_moment = cross(
                core::array::from_fn(|i| center[i] - origin[i]),
                receipt.bearing_impulse,
            );
            near(
                receipt.bearing_angular_impulse[2] - bearing_force_moment[2],
                0.0,
                1e-13,
            );
            near(
                receipt.hydrodynamic.wall_work,
                receipt.hydrodynamic_angular_impulse * receipt.angular_velocity_midpoint,
                1e-14,
            );
            near(
                receipt.kinetic_energy_change,
                0.5 * inertia * (rotor.angular_velocity().powi(2) - spin_before.powi(2)),
                1e-17,
            );
            near(receipt.mechanical_work_residual, 0.0, 1e-14);
        }
    }
}

#[test]
fn swirling_fluid_drives_an_initially_stationary_rotor() {
    let n = 8;
    let origin = center(n, Face3::XMax);
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let r = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let u = cross([0.004, 0.0, 0.0], r);
                grid.set_populations(x, y, z, &equilibrium3(1.0, u));
            }
        }
    }
    let mut rotor =
        WallRotor3::new(&grid, Face3::XMax, origin, 2.6, 50.0, 0.0, n * n, &gate).unwrap();
    let receipt = grid
        .step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &gate)
        .unwrap();
    assert!(rotor.angular_velocity() > 0.0 && rotor.angular_velocity() < 0.004);
    assert!(receipt.hydrodynamic_angular_impulse > 0.0);
    assert!(receipt.kinetic_energy_change > 0.0);
    near(receipt.external_work, 0.0, 0.0);
    near(receipt.mechanical_work_residual, 0.0, 1e-15);
}

#[test]
fn incomplete_or_stale_patches_and_unsupported_domains_are_refused() {
    let n = 8;
    let origin = center(n, Face3::ZMin);
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
    assert!(WallRotor3::new(&grid, Face3::ZMin, origin, 4.0, 30.0, 0.001, n * n, &gate).is_err());
    assert!(
        WallRotor3::new(
            &grid,
            Face3::ZMin,
            origin,
            2.6,
            30.0,
            0.001,
            n * n - 1,
            &gate
        )
        .is_err()
    );
    assert!(WallRotor3::new(&grid, Face3::ZMin, origin, 2.6, 0.0, 0.001, n * n, &gate).is_err());
    let mut rotor =
        WallRotor3::new(&grid, Face3::ZMin, origin, 2.6, 30.0, 0.001, n * n, &gate).unwrap();
    grid.voxelize_sdf(|p| if p == [3.5, 3.5, 1.5] { -1.0 } else { 1.0 });
    let before = snapshot(&grid);
    assert!(
        grid.step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &gate)
            .is_err()
    );
    assert_eq!(snapshot(&grid), before);
    near(rotor.angular_velocity(), 0.001, 0.0);
    near(rotor.angle(), 0.0, 0.0);
    for boundaries in [
        BoundarySpec3::periodic(),
        BoundarySpec3::velocity_pressure_duct_z([0.0; 3], 1.0),
    ] {
        let other = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], boundaries);
        assert!(
            WallRotor3::new(&other, Face3::ZMin, origin, 2.6, 30.0, 0.001, n * n, &gate).is_err()
        );
    }
}

#[test]
fn late_native_stream_refusal_and_cancellation_publish_neither_state() {
    let n = 8;
    let origin = center(n, Face3::ZMin);
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(n, n, n, 1e12, [0.0; 3], closed());
    let mut rotor =
        WallRotor3::new(&grid, Face3::ZMin, origin, 2.6, 30.0, 0.001, n * n, &gate).unwrap();
    // Finite positive-density input survives collision but sends a large
    // negative population to its neighbour, refusing in unpublished streaming.
    let mut bad = equilibrium3(1.0, [0.0; 3]);
    bad[0] += 2.0;
    bad[1] -= 2.0;
    grid.set_populations(3, 3, 3, &bad);
    let before = snapshot(&grid);
    let error = grid
        .step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &gate)
        .unwrap_err();
    assert!(matches!(
        error,
        WallRotorError3::Boundary(BoundaryStepError3::Unphysical { .. })
    ));
    assert_eq!(snapshot(&grid), before);
    near(rotor.angular_velocity(), 0.001, 0.0);
    near(rotor.angle(), 0.0, 0.0);
    grid.set_populations(3, 3, 3, &equilibrium3(1.0, [0.0; 3]));
    let before = snapshot(&grid);
    assert!(
        grid.step_with_finite_inertia_rotor(&mut rotor, 1e9, origin, &gate)
            .is_err()
    );
    assert_eq!(snapshot(&grid), before);
    near(rotor.angular_velocity(), 0.001, 0.0);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        grid.step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &cancelled),
        Err(WallRotorError3::Boundary(BoundaryStepError3::Cancelled))
    ));
    assert_eq!(snapshot(&grid), before);
    near(rotor.angle(), 0.0, 0.0);
    // The next buffer retained after refusal is reusable by the same real step.
    grid.step_with_finite_inertia_rotor(&mut rotor, 0.0, origin, &gate)
        .unwrap();
    assert_ne!(snapshot(&grid), before);
    assert!(rotor.angular_velocity() < 0.001);
}

#[test]
fn one_rotor_batch_is_bitwise_equivalent_to_the_existing_api() {
    let n = 8;
    let gate = CancelGate::new();
    let origin = [-2.0, 1.5, 3.0];
    let center = center(n, Face3::YMax);
    let mut single_grid = BoundaryGrid3::new(n, n, n, 0.8, [1e-6, -2e-6, 3e-6], closed());
    let mut batch_grid = BoundaryGrid3::new(n, n, n, 0.8, [1e-6, -2e-6, 3e-6], closed());
    single_grid.perturb(17, 0.0003);
    batch_grid.perturb(17, 0.0003);
    let mut single_rotor = WallRotor3::new(
        &single_grid,
        Face3::YMax,
        center,
        2.6,
        80.0,
        0.003,
        n * n,
        &gate,
    )
    .unwrap();
    let mut batch_rotors = [single_rotor.clone()];
    for external in [0.0, 0.002, -0.001, 0.0] {
        let single = single_grid
            .step_with_finite_inertia_rotor(&mut single_rotor, external, origin, &gate)
            .unwrap();
        let batch = batch_grid
            .step_with_finite_inertia_rotors(&mut batch_rotors, [external], origin, &gate)
            .unwrap();
        let response = batch.rotors[0];
        assert_eq!(snapshot(&single_grid), snapshot(&batch_grid));
        assert_eq!(
            single_rotor.angular_velocity(),
            batch_rotors[0].angular_velocity()
        );
        assert_eq!(single_rotor.angle(), batch_rotors[0].angle());
        assert_eq!(single.hydrodynamic, response.hydrodynamic);
        assert_eq!(single.stationary_walls, batch.stationary_walls);
        assert_eq!(single.body_force_impulse, batch.body_force_impulse);
        assert_eq!(
            single.body_force_angular_impulse,
            batch.body_force_angular_impulse
        );
        assert_eq!(single.bearing_impulse, response.bearing_impulse);
        assert_eq!(
            single.bearing_angular_impulse,
            response.bearing_angular_impulse
        );
        assert_eq!(
            [
                single.angular_velocity_before,
                single.angular_velocity_midpoint,
                single.angular_velocity_after,
                single.hydrodynamic_angular_impulse,
                single.external_angular_impulse,
                single.kinetic_energy_change,
                single.external_work,
                single.mechanical_work_residual
            ],
            [
                response.angular_velocity_before,
                response.angular_velocity_midpoint,
                response.angular_velocity_after,
                response.hydrodynamic_angular_impulse,
                response.external_angular_impulse,
                response.kinetic_energy_change,
                response.external_work,
                response.mechanical_work_residual
            ]
        );
    }
}

#[test]
fn opposite_rotors_counterrotate_in_one_source_free_fluid_timeline() {
    let n = 8;
    let gate = CancelGate::new();
    let origin = [1.25, 2.75, -0.5];
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
    let inertia = 200.0;
    let mut rotors = [
        WallRotor3::new(
            &grid,
            Face3::ZMax,
            center(n, Face3::ZMax),
            2.6,
            inertia,
            -0.003,
            n * n,
            &gate,
        )
        .unwrap(),
        WallRotor3::new(
            &grid,
            Face3::ZMin,
            center(n, Face3::ZMin),
            2.6,
            inertia,
            0.003,
            n * n,
            &gate,
        )
        .unwrap(),
    ];
    let initial = inventory(&grid, origin);
    let mut housing_linear = [0.0; 3];
    let mut housing_angular = [0.0; 3];
    for _ in 0..16 {
        let before = inventory(&grid, origin);
        let spins_before = rotors.each_ref().map(|rotor| rotor.angular_velocity());
        let exchange = grid
            .step_with_finite_inertia_rotors(&mut rotors, [0.0; 2], origin, &gate)
            .unwrap();
        let after = inventory(&grid, origin);
        near(after.0, before.0, 3e-10);
        near(
            rotors[0].angular_velocity(),
            -rotors[1].angular_velocity(),
            3e-16,
        );
        let spin_impulse = (0..2)
            .map(|i| inertia * (rotors[i].angular_velocity() - spins_before[i]))
            .sum::<f64>();
        for axis in 0..3 {
            let bearing_force: f64 = exchange
                .rotors
                .iter()
                .map(|r| r.bearing_impulse[axis])
                .sum();
            let bearing_moment: f64 = exchange
                .rotors
                .iter()
                .map(|r| r.bearing_angular_impulse[axis])
                .sum();
            let force = exchange.stationary_walls.wall_impulse[axis] - bearing_force;
            let torque = exchange.stationary_walls.wall_angular_impulse[axis] - bearing_moment;
            housing_linear[axis] += force;
            housing_angular[axis] += torque;
            near(after.1[axis] - before.1[axis] + force, 0.0, 2e-12);
            near(
                after.2[axis] - before.2[axis]
                    + torque
                    + if axis == 2 { spin_impulse } else { 0.0 },
                0.0,
                8e-12,
            );
        }
        assert_eq!(exchange.body_force_impulse, [0.0; 3]);
        assert_eq!(exchange.body_force_angular_impulse, [0.0; 3]);
        for response in exchange.rotors {
            near(response.external_work, 0.0, 0.0);
            near(response.mechanical_work_residual, 0.0, 2e-14);
        }
    }
    let final_state = inventory(&grid, origin);
    for axis in 0..3 {
        near(
            final_state.1[axis] - initial.1[axis] + housing_linear[axis],
            0.0,
            2e-11,
        );
        near(
            final_state.2[axis] - initial.2[axis]
                + housing_angular[axis]
                + if axis == 2 {
                    inertia * (rotors[0].angular_velocity() + rotors[1].angular_velocity())
                } else {
                    0.0
                },
            0.0,
            8e-11,
        );
    }
    assert!(
        rotors
            .iter()
            .all(|rotor| rotor.angular_velocity().abs() < 0.003)
    );
}

#[test]
fn rotors_on_three_spatial_axes_share_one_housing_and_force_exchange() {
    let n = 8;
    let gate = CancelGate::new();
    let origin = [-1.1, 2.3, 9.7];
    let force = [1e-6, 2e-6, -1e-6];
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, force, closed());
    grid.perturb(81, 0.0005);
    let faces = [Face3::XMin, Face3::YMax, Face3::ZMin];
    let inertias = [53.0, 39.0, 67.0];
    let spins = [0.003, -0.002, 0.001];
    let external = [0.004, -0.003, 0.005];
    let mut rotors = core::array::from_fn(|i| {
        WallRotor3::new(
            &grid,
            faces[i],
            center(n, faces[i]),
            1.6,
            inertias[i],
            spins[i],
            n * n,
            &gate,
        )
        .unwrap()
    });
    for _ in 0..4 {
        let before = inventory(&grid, origin);
        let spins_before = rotors.each_ref().map(|rotor| rotor.angular_velocity());
        let exchange = grid
            .step_with_finite_inertia_rotors(&mut rotors, external, origin, &gate)
            .unwrap();
        let after = inventory(&grid, origin);
        near(after.0, before.0, 3e-10);
        for axis in 0..3 {
            let rotor_impulse =
                inertias[axis] * (rotors[axis].angular_velocity() - spins_before[axis]);
            let bearing_force: f64 = exchange
                .rotors
                .iter()
                .map(|r| r.bearing_impulse[axis])
                .sum();
            let bearing_moment: f64 = exchange
                .rotors
                .iter()
                .map(|r| r.bearing_angular_impulse[axis])
                .sum();
            near(
                exchange.body_force_impulse[axis],
                (n * n * n) as f64 * force[axis],
                3e-17,
            );
            near(
                after.1[axis] - before.1[axis] + exchange.stationary_walls.wall_impulse[axis]
                    - bearing_force,
                exchange.body_force_impulse[axis],
                3e-12,
            );
            near(
                after.2[axis] - before.2[axis]
                    + rotor_impulse
                    + exchange.stationary_walls.wall_angular_impulse[axis]
                    - bearing_moment,
                exchange.body_force_angular_impulse[axis] + external[axis],
                1e-11,
            );
            near(
                rotor_impulse,
                external[axis] + exchange.rotors[axis].hydrodynamic_angular_impulse,
                2e-12,
            );
            near(exchange.rotors[axis].mechanical_work_residual, 0.0, 2e-14);
        }
        let measured = exchange.stationary_walls.measured_links
            + exchange
                .rotors
                .iter()
                .map(|r| r.hydrodynamic.measured_links)
                .sum::<usize>();
        assert_eq!(measured, grid.boundary_links().len());
    }
}

#[test]
fn overlapping_patches_or_one_late_rotor_refusal_leave_every_state_intact() {
    let n = 8;
    let gate = CancelGate::new();
    let origin = [0.0; 3];
    let mut grid = BoundaryGrid3::new(n, n, n, 0.8, [0.0; 3], closed());
    let left = WallRotor3::new(
        &grid,
        Face3::ZMin,
        [2.5, 4.0, 0.0],
        1.2,
        30.0,
        0.003,
        n * n,
        &gate,
    )
    .unwrap();
    let right = WallRotor3::new(
        &grid,
        Face3::ZMin,
        [5.5, 4.0, 0.0],
        1.2,
        40.0,
        -0.002,
        n * n,
        &gate,
    )
    .unwrap();
    let mut overlap = [left.clone(), left.clone()];
    let before = snapshot(&grid);
    let error = grid
        .step_with_finite_inertia_rotors(&mut overlap, [0.0; 2], origin, &gate)
        .unwrap_err();
    assert!(matches!(
        error,
        WallRotorError3::Invalid("rotor patches overlap")
    ));
    assert_eq!(snapshot(&grid), before);
    assert!(
        overlap
            .iter()
            .all(|rotor| rotor.angle() == 0.0 && rotor.angular_velocity() == 0.003)
    );

    let mut disjoint = [left, right];
    let error = grid
        .step_with_finite_inertia_rotors(&mut disjoint, [0.0, 1e9], origin, &gate)
        .unwrap_err();
    assert!(matches!(error, WallRotorError3::Invalid(_)));
    assert_eq!(snapshot(&grid), before);
    assert_eq!(
        disjoint.each_ref().map(|rotor| rotor.angular_velocity()),
        [0.003, -0.002]
    );
    assert_eq!(disjoint.each_ref().map(|rotor| rotor.angle()), [0.0; 2]);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        grid.step_with_finite_inertia_rotors(&mut disjoint, [0.0; 2], origin, &cancelled),
        Err(WallRotorError3::Boundary(BoundaryStepError3::Cancelled))
    ));
    assert_eq!(snapshot(&grid), before);
    let mut empty: [WallRotor3; 0] = [];
    assert!(
        grid.step_with_finite_inertia_rotors(&mut empty, [], origin, &gate)
            .is_err()
    );
    assert_eq!(snapshot(&grid), before);

    // Two separately driven patches on the same face are a valid joint step.
    let exchange = grid
        .step_with_finite_inertia_rotors(&mut disjoint, [0.0; 2], origin, &gate)
        .unwrap();
    assert_ne!(snapshot(&grid), before);
    assert!(disjoint[0].angular_velocity() < 0.003);
    assert!(disjoint[1].angular_velocity() > -0.002);
    assert_eq!(
        exchange.stationary_walls.measured_links
            + exchange
                .rotors
                .iter()
                .map(|r| r.hydrodynamic.measured_links)
                .sum::<usize>(),
        grid.boundary_links().len()
    );
}
