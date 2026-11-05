//! G0/G1/G3/G4/G5: actual-density body acceleration on the dense native owner.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{BoundaryStepError3, CollisionModel3, E3, Q3};
use fs_lbm::{BoundaryGrid3, BoundarySpec3, equilibrium3};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn cells(grid: &BoundaryGrid3) -> impl Iterator<Item = [usize; 3]> + '_ {
    let [nx, ny, nz] = grid.dimensions();
    (0..nz)
        .flat_map(move |z| (0..ny).flat_map(move |y| (0..nx).map(move |x| [x, y, z])))
        .filter(move |&[x, y, z]| !grid.is_solid(x, y, z))
}

fn moments(grid: &BoundaryGrid3, [x, y, z]: [usize; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let populations = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let density = populations.iter().sum::<f64>();
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let momentum: [f64; 3] = core::array::from_fn(|a| {
        populations
            .iter()
            .zip(directions)
            .map(|(f, e)| f * e[a])
            .sum::<f64>()
            + 0.5 * force[a]
    });
    let velocity = momentum.map(|p| p / density);
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            populations
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                - density * velocity[a] * velocity[b]
                - if a == b { density / 3.0 } else { 0.0 }
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
        })
    });
    (density, momentum, stress)
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    cells(grid)
        .map(|[x, y, z]| {
            (
                grid.populations(x, y, z).map(f64::to_bits),
                grid.force_density(x, y, z).map(f64::to_bits),
            )
        })
        .collect()
}

fn models() -> [CollisionModel3; 4] {
    [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 0.9,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 0.9,
            fourth_order_rate: 1.1,
        },
    ]
}

#[test]
fn acceleration_and_jerk_response_do_not_depend_on_density_or_collision_rate() {
    let gate = CancelGate::new_clock_free();
    for model in models() {
        for density in [0.4, 1.0, 2.3] {
            let initial = [0.011, -0.009, 0.006];
            let mut grid = BoundaryGrid3::with_collision_model(
                4,
                4,
                4,
                model,
                [0.0; 3],
                BoundarySpec3::periodic(),
            );
            grid.initialize_uniform_equilibrium(density, initial, &gate)
                .unwrap();
            let mut previous = [0.0007, -0.0012, 0.0009];
            grid.rebind_body_acceleration(previous, &gate).unwrap();
            let mut expected = initial;
            for acceleration in [previous, [-0.0003, 0.0004, 0.0008], [0.0; 3], [0.0; 3]] {
                grid.step_with_endpoint_body_acceleration(acceleration, &gate)
                    .unwrap();
                let (rho, momentum, _) = moments(&grid, [1, 2, 3]);
                let force = grid.force_density(1, 2, 3);
                for a in 0..3 {
                    expected[a] += 0.5 * (previous[a] + acceleration[a]);
                    close(momentum[a] / rho, expected[a], 8e-16);
                    assert_eq!(force[a], rho * acceleration[a]);
                }
                previous = acceleration;
            }
        }
    }
}

#[test]
fn rebinding_replaces_previous_sources_without_erasing_physical_stress_or_impulse() {
    let gate = CancelGate::new_clock_free();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    let mut populations = equilibrium3(1.6, [0.03, -0.02, 0.01]);
    populations[1] += 0.002;
    populations[2] += 0.002;
    populations[0] -= 0.004;
    populations[7] += 0.001;
    populations[8] += 0.001;
    populations[9] -= 0.001;
    populations[10] -= 0.001;
    for cell in cells(&grid).collect::<Vec<_>>() {
        grid.set_populations(cell[0], cell[1], cell[2], &populations);
    }
    grid.rebind_uniform_force([0.0002, -0.0003, 0.0001], &gate)
        .unwrap();
    grid.rebind_force_field(|[x, _, _]| [x as f64 * 0.0001, 0.0, 0.0]);
    let before: Vec<_> = cells(&grid).map(|cell| moments(&grid, cell)).collect();
    for acceleration in [[0.001, -0.002, 0.003], [-0.002, 0.001, -0.001], [0.0; 3]] {
        grid.rebind_body_acceleration(acceleration, &gate).unwrap();
        assert_eq!(grid.uniform_force_density(), [0.0; 3]);
        for (cell, original) in cells(&grid).zip(&before) {
            let (rho, momentum, stress) = moments(&grid, cell);
            close(rho, original.0, 2e-15);
            for a in 0..3 {
                close(momentum[a], original.1[a], 2e-16);
                close(
                    grid.force_density(cell[0], cell[1], cell[2])[a],
                    rho * acceleration[a],
                    1e-17,
                );
                for b in 0..3 {
                    close(stress[a][b], original.2[a][b], 3e-16);
                }
            }
        }
    }
}

#[test]
fn transported_density_enters_the_endpoint_source_without_a_fictitious_reframe() {
    let gate = CancelGate::new_clock_free();
    let gradient = [0.008, -0.006, 0.004];
    let start = [0.004, -0.003, 0.002];
    let end = [-0.002, 0.001, 0.005];
    for model in models() {
        let mut grid = BoundaryGrid3::with_collision_model(
            8,
            8,
            8,
            model,
            [0.0; 3],
            BoundarySpec3::periodic(),
        );
        for cell in cells(&grid).collect::<Vec<_>>() {
            let density = 1.3
                + (0..3)
                    .map(|a| gradient[a] * (cell[a] as f64 - 4.0))
                    .sum::<f64>();
            grid.set_populations(cell[0], cell[1], cell[2], &equilibrium3(density, [0.0; 3]));
        }
        grid.rebind_body_acceleration(start, &gate).unwrap();
        grid.step_with_endpoint_body_acceleration(end, &gate)
            .unwrap();
        // The independent D3Q19 sums on an affine density field at rest give
        // rho1 = rho0 - (grad rho).a0/2 and raw j1 = -grad rho/3 + rho0*a0/2.
        // A held rho0*a1 or a post-stream reframe changes the physical impulse.
        let expected_density = 1.3 - 0.5 * (0..3).map(|a| gradient[a] * start[a]).sum::<f64>();
        let (density, momentum, _) = moments(&grid, [4, 4, 4]);
        close(density, expected_density, 2e-15);
        for a in 0..3 {
            let expected = -gradient[a] / 3.0 + 0.5 * (1.3 * start[a] + expected_density * end[a]);
            close(momentum[a], expected, 3e-16);
            assert_eq!(grid.force_density(4, 4, 4)[a], density * end[a]);
        }
    }
}

#[test]
fn smooth_acceleration_history_has_second_order_physical_velocity_error() {
    let gate = CancelGate::new_clock_free();
    let omega: f64 = 2.1;
    let amplitude = 0.1;
    let exact = amplitude * omega.sin() / omega;
    let errors = [10, 20, 40].map(|steps| {
        let dt = 1.0 / f64::from(steps);
        let acceleration = |time: f64| [amplitude * (omega * time).cos() * dt * dt, 0.0, 0.0];
        let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
        grid.rebind_body_acceleration(acceleration(0.0), &gate)
            .unwrap();
        for step in 1..=steps {
            grid.step_with_endpoint_body_acceleration(acceleration(f64::from(step) * dt), &gate)
                .unwrap();
        }
        let (rho, momentum, _) = moments(&grid, [0, 0, 0]);
        (momentum[0] / rho / dt - exact).abs()
    });
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.98 && pair[0] / pair[1] < 4.02,
            "{errors:?}"
        );
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn totals(grid: &BoundaryGrid3, origin: [f64; 3]) -> [[f64; 3]; 4] {
    let mut sums = [[0.0; 3]; 4];
    for cell in cells(grid) {
        let momentum = moments(grid, cell).1;
        let force = grid.force_density(cell[0], cell[1], cell[2]);
        let r = core::array::from_fn(|a| cell[a] as f64 + 0.5 - origin[a]);
        for (sum, value) in
            sums.iter_mut()
                .zip([momentum, cross(r, momentum), force, cross(r, force)])
        {
            for a in 0..3 {
                sum[a] += value[a];
            }
        }
    }
    sums
}

#[test]
fn periodic_transport_retains_trapezoidal_total_impulse_and_each_local_density() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xACCE0);
    let mut serial = BoundaryGrid3::new(8, 4, 4, 0.85, [0.0; 3], BoundarySpec3::periodic());
    for cell in cells(&serial).collect::<Vec<_>>() {
        let angle = std::f64::consts::TAU * cell[0] as f64 / 8.0;
        serial.set_populations(
            cell[0],
            cell[1],
            cell[2],
            &equilibrium3(
                0.8 + 0.06 * angle.cos(),
                [0.012 + 0.004 * angle.sin(), -0.003, 0.002],
            ),
        );
    }
    serial
        .rebind_body_acceleration([0.001, -0.002, 0.0015], &gate)
        .unwrap();
    let mut pooled = serial.clone();
    for acceleration in [[-0.003, 0.002, 0.001], [0.001, -0.002, -0.003], [0.0; 3]] {
        let before = totals(&serial, [0.0; 3]);
        serial
            .step_with_endpoint_body_acceleration(acceleration, &gate)
            .unwrap();
        pooled
            .step_pooled_with_endpoint_body_acceleration(&pool, &gate, acceleration)
            .unwrap();
        assert_eq!(snapshot(&serial), snapshot(&pooled));
        let after = totals(&serial, [0.0; 3]);
        for a in 0..3 {
            close(
                after[0][a] - before[0][a],
                0.5 * (before[2][a] + after[2][a]),
                1e-12,
            );
        }
        for cell in cells(&serial) {
            let density = serial
                .populations(cell[0], cell[1], cell[2])
                .iter()
                .sum::<f64>();
            assert_eq!(
                serial.force_density(cell[0], cell[1], cell[2]),
                acceleration.map(|a| density * a)
            );
        }
    }
}

#[test]
fn density_source_and_native_wall_loads_close_linear_and_angular_impulse_with_pooled_bits() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(3, 0xACCE1);
    let origin = [-0.3, 0.7, 1.2];
    let mut serial =
        BoundaryGrid3::new(8, 8, 8, 0.9, [0.0; 3], BoundarySpec3::lid_cavity([0.0; 3]));
    serial.voxelize_sdf(|p| {
        (p[0] - 3.0)
            .abs()
            .max((p[1] - 4.0).abs())
            .max((p[2] - 3.0).abs())
            - 0.8
    });
    for cell in cells(&serial).collect::<Vec<_>>() {
        let rho = 0.9 + 0.01 * cell[0] as f64 - 0.005 * cell[1] as f64;
        serial.set_populations(
            cell[0],
            cell[1],
            cell[2],
            &equilibrium3(rho, [0.003, -0.002, 0.001]),
        );
    }
    serial
        .rebind_body_acceleration([0.0007, -0.0012, 0.0009], &gate)
        .unwrap();
    let mut pooled = serial.clone();
    let links = serial.boundary_links();
    let mass = serial.total_mass();
    for acceleration in [
        [0.0003, -0.0002, 0.0001],
        [-0.0005, 0.0004, -0.0006],
        [0.0; 3],
    ] {
        let before = totals(&serial, origin);
        let exchange = serial
            .step_with_endpoint_body_acceleration_and_wall_momentum(
                acceleration,
                &links,
                origin,
                &gate,
            )
            .unwrap();
        let parallel = pooled
            .step_pooled_with_endpoint_body_acceleration_and_wall_momentum(
                &pool,
                &gate,
                acceleration,
                &links,
                origin,
            )
            .unwrap();
        let after = totals(&serial, origin);
        assert_eq!(exchange, parallel);
        assert_eq!(snapshot(&serial), snapshot(&pooled));
        close(serial.total_mass(), mass, 2e-12);
        for a in 0..3 {
            close(
                after[0][a] - before[0][a] + exchange.wall_impulse[a],
                0.5 * (before[2][a] + after[2][a]),
                2e-12,
            );
            close(
                after[1][a] - before[1][a] + exchange.wall_angular_impulse[a],
                0.5 * (before[3][a] + after[3][a]),
                5e-12,
            );
        }
    }
}

#[test]
fn late_endpoint_refusal_and_cancellation_preserve_accepted_state() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xACCE2);
    for pooled in [false, true] {
        let advance = |grid: &mut BoundaryGrid3, acceleration, gate: &CancelGate| {
            if pooled {
                grid.step_pooled_with_endpoint_body_acceleration(&pool, gate, acceleration)
            } else {
                grid.step_with_endpoint_body_acceleration(acceleration, gate)
            }
        };
        let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
        grid.initialize_uniform_equilibrium(2.3, [0.01, 0.0, 0.0], &gate)
            .unwrap();
        // No field has been allocated/bound. Failed first-use publication
        // must retain the original zero source as well as population bits.
        let before = snapshot(&grid);
        assert!(matches!(
            advance(&mut grid, [f64::MAX, 0.0, 0.0], &gate),
            Err(BoundaryStepError3::Unphysical { .. })
        ));
        assert_eq!(snapshot(&grid), before);
        assert!(advance(&mut grid, [f64::NAN, 0.0, 0.0], &gate).is_err());
        assert_eq!(snapshot(&grid), before);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert_eq!(
            advance(&mut grid, [0.001, 0.0, 0.0], &cancelled),
            Err(BoundaryStepError3::Cancelled)
        );
        assert_eq!(snapshot(&grid), before);
        // A refused proposal must not lock topology as if a step succeeded.
        grid.voxelize_sdf(|_| 1.0);
        advance(&mut grid, [0.001, 0.0, 0.0], &gate).unwrap();
        assert_ne!(snapshot(&grid), before);
    }
}
