//! Physical Guo moments under held uniform loads, including retained stress.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, E3, Q3, equilibrium3};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn populations(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut out = Vec::new();
    let [nx, ny, nz] = grid.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                out.push(grid.populations(x, y, z));
            }
        }
    }
    out
}
fn moments(grid: &BoundaryGrid3, cell: [usize; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let [x, y, z] = cell;
    let f = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
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
fn total_momentum(grid: &BoundaryGrid3) -> [f64; 3] {
    let [nx, ny, nz] = grid.dimensions();
    let mut out = [0.0; 3];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let p = moments(grid, [x, y, z]).1;
                for a in 0..3 {
                    out[a] += p[a];
                }
            }
        }
    }
    out
}

#[test]
fn uniform_rebind_preserves_local_momentum_stress_and_existing_spatial_force() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    let mut f = equilibrium3(1.3, [0.03, -0.02, 0.01]);
    // Retained nonequilibrium normal stress distinguishes a source-frame change
    // from resetting the evolved fluid to equilibrium.
    f[1] += 0.002;
    f[2] += 0.002;
    f[0] -= 0.004;
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                grid.set_populations(x, y, z, &f);
            }
        }
    }
    let spatial = |[x, y, z]: [usize; 3]| [x as f64 * 1e-5, y as f64 * -2e-5, z as f64 * 3e-5];
    grid.rebind_force_field(spatial);
    let original = moments(&grid, [2, 1, 3]);
    for force in [[0.001, -0.002, 0.003], [-0.004, 0.001, -0.002], [0.0; 3]] {
        grid.rebind_uniform_force(force, &gate).unwrap();
        let next = moments(&grid, [2, 1, 3]);
        near(next.0, original.0, 2e-15);
        for a in 0..3 {
            near(next.1[a], original.1[a], 2e-15);
            for b in 0..3 {
                near(next.2[a][b], original.2[a][b], 2e-15);
            }
            near(
                grid.force_density(2, 1, 3)[a],
                force[a] + spatial([2, 1, 3])[a],
                0.0,
            );
        }
        let accepted = populations(&grid);
        grid.rebind_uniform_force(force, &gate).unwrap();
        assert_eq!(populations(&grid), accepted);
    }
}

#[test]
fn physical_initializer_and_switched_load_apply_only_elapsed_impulse() {
    let gate = CancelGate::new();
    let first = [0.0002, -0.0001, 0.0003];
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, first, BoundarySpec3::periodic());
    let velocity = [0.02, -0.01, 0.03];
    grid.initialize_uniform_equilibrium(1.7, velocity, &gate)
        .unwrap();
    let initial = moments(&grid, [1, 2, 3]);
    near(initial.0, 1.7, 2e-15);
    for a in 0..3 {
        near(initial.1[a], 1.7 * velocity[a], 2e-15);
        for b in 0..3 {
            near(
                initial.2[a][b],
                1.7 * velocity[a] * velocity[b] + if a == b { 1.7 / 3.0 } else { 0.0 },
                2e-15,
            );
        }
    }
    for force in [first, [-0.0003, 0.0002, -0.0001], [0.0; 3]] {
        let before = total_momentum(&grid);
        grid.rebind_uniform_force(force, &gate).unwrap();
        let rebound = total_momentum(&grid);
        for a in 0..3 {
            near(rebound[a], before[a], 2e-13);
        }
        grid.step();
        let after = total_momentum(&grid);
        for a in 0..3 {
            near(after[a] - before[a], 64.0 * force[a], 2e-13);
        }
    }
}

#[test]
fn late_bad_cell_and_cancel_refuse_while_open_force_rebinding_preserves_momentum() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    grid.initialize_uniform_equilibrium(1.0, [0.02, -0.01, 0.03], &gate)
        .unwrap();
    grid.set_populations(3, 3, 3, &[-1.0; Q3]);
    let before = populations(&grid);
    assert_eq!(
        grid.rebind_uniform_force([0.001; 3], &gate),
        Err(BoundaryForcingError3::Unphysical { cell: [3, 3, 3] })
    );
    assert_eq!(populations(&grid), before);
    assert_eq!(grid.force_density(0, 0, 0), [0.0; 3]);
    gate.request();
    assert_eq!(
        grid.initialize_uniform_equilibrium(1.0, [0.0; 3], &gate),
        Err(BoundaryForcingError3::Cancelled)
    );
    assert_eq!(
        grid.rebind_uniform_force([0.001; 3], &gate),
        Err(BoundaryForcingError3::Cancelled)
    );
    assert_eq!(populations(&grid), before);
    let mut open = BoundaryGrid3::new(
        4,
        4,
        4,
        0.9,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(1.0, 1.0),
    );
    let before: Vec<_> = (0..4).map(|z| moments(&open, [1, 1, z])).collect();
    open.rebind_uniform_force([0.0, 0.0, 0.001], &CancelGate::new())
        .unwrap();
    for (z, expected) in before.iter().enumerate() {
        let actual = moments(&open, [1, 1, z]);
        near(actual.0, expected.0, 2e-15);
        for a in 0..3 {
            near(actual.1[a], expected.1[a], 2e-15);
            for b in 0..3 {
                near(actual.2[a][b], expected.2[a][b], 2e-15);
            }
        }
    }
}
