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
