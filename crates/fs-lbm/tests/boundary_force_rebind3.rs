//! Spatial held-force admission must preserve physical state until real time elapses.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, CollisionModel3, E3, Q3};
use fs_lbm::{BoundaryGrid3, BoundarySpec3, CS2, equilibrium3};

fn cells(grid: &BoundaryGrid3) -> Vec<[usize; 3]> {
    let [nx, ny, nz] = grid.dimensions();
    (0..nz)
        .flat_map(|z| (0..ny).flat_map(move |y| (0..nx).map(move |x| [x, y, z])))
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn field([x, y, z]: [usize; 3]) -> [f64; 3] {
    [
        1e-5 * (1 + x + 2 * z) as f64,
        -2e-5 * (1 + y) as f64,
        3e-5 * (1 + z + y) as f64,
    ]
}

// Independent physical moments: do not use the production velocity or source transform.
fn moments(grid: &BoundaryGrid3, [x, y, z]: [usize; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let f = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let rho = f.iter().sum::<f64>();
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let p = core::array::from_fn(|a| {
        f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]
    });
    let u = p.map(|p| p / rho);
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                - rho * u[a] * u[b]
                - if a == b { rho * CS2 } else { 0.0 }
                + 0.5 * (u[a] * force[b] + force[a] * u[b])
        })
    });
    (rho, p, stress)
}

fn momentum(grid: &BoundaryGrid3) -> [f64; 3] {
    cells(grid).into_iter().fold([0.0; 3], |sum, cell| {
        let p = moments(grid, cell).1;
        core::array::from_fn(|a| sum[a] + p[a])
    })
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    cells(grid)
        .into_iter()
        .filter(|&[x, y, z]| !grid.is_solid(x, y, z))
        .map(|[x, y, z]| {
            (
                grid.populations(x, y, z).map(f64::to_bits),
                grid.force_density(x, y, z).map(f64::to_bits),
            )
        })
        .collect()
}

#[test]
fn spatial_rebind_preserves_full_physical_stress_and_x_fastest_sampling() {
    let gate = CancelGate::new_clock_free();
    let uniform = [2e-4, -1e-4, 3e-4];
    let mut grid = BoundaryGrid3::new(8, 4, 8, 0.9, uniform, BoundarySpec3::periodic());
    grid.voxelize_sdf(|p| (p[0] - 3.5).abs() - 0.25);
    grid.set_force_field(field);
    for [x, y, z] in cells(&grid) {
        if grid.is_solid(x, y, z) {
            continue;
        }
        let mut f = equilibrium3(
            1.3 + y as f64 * 0.07,
            [0.01, -0.02, 0.005 + z as f64 * 0.001],
        );
        // Genuine normal and shear nonequilibrium stress at zero added mass/momentum.
        f[1] += 0.002;
        f[2] += 0.002;
        f[0] -= 0.004;
        f[7] += 0.001;
        f[8] += 0.001;
        f[9] -= 0.001;
        f[10] -= 0.001;
        grid.set_populations(x, y, z, &f);
    }
    let all = cells(&grid);
    let fluid: Vec<_> = all
        .iter()
        .copied()
        .filter(|&[x, y, z]| !grid.is_solid(x, y, z))
        .collect();
    let before: Vec<_> = fluid.iter().map(|&cell| moments(&grid, cell)).collect();
    let mut sampled = Vec::new();
    grid.rebind_force_field_with_gate(
        |cell| {
            sampled.push(cell);
            field(cell).map(|f| -0.7 * f)
        },
        &gate,
    )
    .unwrap();
    assert_eq!(sampled, fluid);
    assert_eq!(grid.uniform_force_density(), uniform);
    for (&cell, expected) in fluid.iter().zip(before) {
        let actual = moments(&grid, cell);
        close(actual.0, expected.0, 4e-15);
        for a in 0..3 {
            close(actual.1[a], expected.1[a], 4e-15);
            for b in 0..3 {
                close(actual.2[a][b], expected.2[a][b], 4e-15);
            }
        }
    }
    let rebound = snapshot(&grid);
    grid.rebind_force_field_with_gate(|cell| field(cell).map(|f| -0.7 * f), &gate)
        .unwrap();
    assert_eq!(
        snapshot(&grid),
        rebound,
        "unchanged force must preserve population bits"
    );
    grid.set_force_field_scale(0.5);
    for [x, y, z] in fluid {
        for a in 0..3 {
            close(
                grid.force_density(x, y, z)[a],
                uniform[a] - 0.35 * field([x, y, z])[a],
                2e-19,
            );
        }
    }
}

#[test]
fn switched_spatial_loads_supply_only_elapsed_collision_impulses() {
    let gate = CancelGate::new_clock_free();
    let uniform = [2e-5, -1e-5, 3e-5];
    for collision in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 0.9,
        },
    ] {
        let mut grid = BoundaryGrid3::with_collision_model(
            8,
            4,
            4,
            collision,
            [0.0; 3],
            BoundarySpec3::periodic(),
        );
        grid.initialize_uniform_equilibrium(1.4, [0.012, -0.004, 0.003], &gate)
            .unwrap();
        grid.rebind_uniform_force(uniform, &gate).unwrap();
        let mass = grid.total_mass();
        let mut expected = momentum(&grid);
        for scale in [1.0, -0.7, 0.0] {
            grid.rebind_force_field_with_gate(|cell| field(cell).map(|f| scale * f), &gate)
                .unwrap();
            for a in 0..3 {
                close(momentum(&grid)[a], expected[a], 3e-13);
            }
            let impulse = cells(&grid).into_iter().fold([0.0; 3], |sum, cell| {
                core::array::from_fn(|a| sum[a] + uniform[a] + scale * field(cell)[a])
            });
            for _ in 0..2 {
                grid.step();
                for a in 0..3 {
                    expected[a] += impulse[a];
                    close(momentum(&grid)[a], expected[a], 8e-13);
                }
                close(grid.total_mass(), mass, 64.0 * f64::EPSILON * mass.abs());
            }
        }
    }
}

#[test]
fn typed_late_refusal_and_cancellation_keep_accepted_state_and_retry_intact() {
    let gate = CancelGate::new_clock_free();
    let make = || {
        let mut grid = BoundaryGrid3::new(8, 4, 4, 0.8, [0.0; 3], BoundarySpec3::periodic());
        grid.rebind_force_field_with_gate(field, &gate).unwrap();
        grid.step();
        grid
    };
    let mut grid = make();
    let mut control = make();
    let before = snapshot(&grid);
    assert_eq!(
        grid.rebind_force_field_with_gate(
            |cell| {
                if cell == [7, 3, 3] {
                    [f64::NAN, 0.0, 0.0]
                } else {
                    field(cell).map(|f| -f)
                }
            },
            &gate
        ),
        Err(BoundaryForcingError3::Unphysical { cell: [7, 3, 3] })
    );
    assert_eq!(snapshot(&grid), before);
    let cancelled = CancelGate::new_clock_free();
    assert_eq!(
        grid.rebind_force_field_with_gate(
            |cell| {
                if cell == [4, 0, 0] {
                    cancelled.request();
                }
                field(cell)
            },
            &cancelled
        ),
        Err(BoundaryForcingError3::Cancelled)
    );
    assert_eq!(snapshot(&grid), before);
    assert_eq!(
        grid.rebind_force_field_with_gate(
            |_| panic!("cancelled callback must not run"),
            &cancelled
        ),
        Err(BoundaryForcingError3::Cancelled)
    );
    let last = grid.populations(7, 3, 3);
    grid.set_populations(7, 3, 3, &[0.0; Q3]);
    let invalid = snapshot(&grid);
    assert_eq!(
        grid.rebind_force_field_with_gate(|cell| field(cell).map(|f| -f), &gate),
        Err(BoundaryForcingError3::Unphysical { cell: [7, 3, 3] })
    );
    assert_eq!(snapshot(&grid), invalid);
    grid.set_populations(7, 3, 3, &last);
    for grid in [&mut grid, &mut control] {
        grid.rebind_force_field_with_gate(|cell| field(cell).map(|f| -f), &gate)
            .unwrap();
        grid.step();
    }
    assert_eq!(
        snapshot(&grid),
        snapshot(&control),
        "refused scratch must not affect an accepted retry"
    );
}
