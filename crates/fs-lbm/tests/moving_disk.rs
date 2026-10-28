//! G0/G1/G3: moving analytic geometry and native population/body momentum.

use fs_lbm::core2::moving_disk::{MovingDisk2, disk_cell_fraction};
use fs_lbm::core2::{Cell, Grid, PartialSaturationBody2};

const E: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn body(center: [f64; 2], velocity: [f64; 2], angular_velocity: f64) -> PartialSaturationBody2 {
    PartialSaturationBody2 {
        mass: 7.0,
        inertia: 5.0,
        center,
        velocity,
        angular_velocity,
    }
}

fn totals(grid: &Grid, body: PartialSaturationBody2) -> [f64; 4] {
    let mut total = [
        body.mass,
        body.mass * body.velocity[0],
        body.mass * body.velocity[1],
        body.inertia * body.angular_velocity + body.center[0] * body.mass * body.velocity[1]
            - body.center[1] * body.mass * body.velocity[0],
    ];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            if grid.flags[i] != Cell::Fluid {
                continue;
            }
            for (population, direction) in grid.f[i].iter().zip(E) {
                let px = population * direction[0];
                let py = population * direction[1];
                total[0] += population;
                total[1] += px;
                total[2] += py;
                total[3] += x as f64 * py - y as f64 * px;
            }
        }
    }
    total
}

#[test]
fn analytic_raster_resolves_small_disks_tangencies_and_translated_area() {
    let radius = 0.25;
    let area = std::f64::consts::PI * radius * radius;
    close(disk_cell_fraction(radius, [0.0, 0.0]), area, 2e-16);
    close(disk_cell_fraction(radius, [0.5, 0.5]), 0.25 * area, 2e-16);
    close(disk_cell_fraction(radius, [0.5, 0.0]), 0.5 * area, 2e-16);
    close(
        disk_cell_fraction(0.5, [0.0, 0.0]),
        0.25 * std::f64::consts::PI,
        2e-16,
    );
    assert_eq!(disk_cell_fraction(0.5, [1.0, 0.0]), 0.0);
    assert_eq!(disk_cell_fraction(0.75, [0.0, 0.0]), 1.0);
    assert_eq!(disk_cell_fraction(0.5, [2.0, 2.0]), 0.0);
    for radius in [0.25, 0.5, 1.0, 2.35, 6.7] {
        for shift in 0..24 {
            let center = [12.0 + 0.037 * shift as f64, 11.0 - 0.023 * shift as f64];
            let mut sum = 0.0;
            for y in 0..26 {
                for x in 0..26 {
                    let offset = [x as f64 - center[0], y as f64 - center[1]];
                    let fraction = disk_cell_fraction(radius, offset);
                    assert!((0.0..=1.0).contains(&fraction));
                    close(
                        fraction,
                        disk_cell_fraction(radius, [offset[1], -offset[0]]),
                        2e-13,
                    );
                    sum += fraction;
                }
            }
            close(sum, std::f64::consts::PI * radius * radius, 5e-12);
        }
    }
}

#[test]
fn common_translation_moves_the_actual_geometry_without_ghost_coupling() {
    let mut grid = Grid::uniform(40, 32, 0.8);
    let velocity = [0.04, 0.015];
    grid.f
        .fill(fs_lbm::equilibrium(1.0, velocity[0], velocity[1]));
    let initial = body([10.0, 11.0], velocity, 0.0);
    let mut disk = MovingDisk2::new(&grid, 2.1, initial, 0.3);
    assert_eq!(disk.solid_fraction(10, 11), 1.0);
    let initial_mass = grid.total_mass();
    for _ in 0..180 {
        let receipt = disk.step(&mut grid, [0.0; 2], 0.0);
        for residual in receipt.exchange.momentum_residual {
            close(residual, 0.0, 2e-13);
        }
    }
    for axis in 0..2 {
        close(disk.body().velocity[axis], velocity[axis], 1e-12);
        close(
            disk.body().center[axis],
            initial.center[axis] + 180.0 * velocity[axis],
            2e-11,
        );
    }
    close(disk.angle(), 0.3, 1e-12);
    close(grid.total_mass(), initial_mass, 2e-9);
    assert_eq!(disk.solid_fraction(10, 11), 0.0);
    assert_eq!(disk.solid_fraction(17, 14), 1.0);
    // A disturbance at the old occupied location must not exert an immediate
    // impulse on the translated body through an uncleared fraction cache.
    let old = grid.idx(10, 11);
    grid.f[old] = fs_lbm::equilibrium(1.0, velocity[0] + 0.01, velocity[1]);
    let before = disk.body();
    let receipt = disk.step(&mut grid, [0.0; 2], 0.0);
    for axis in 0..2 {
        close(disk.body().velocity[axis], before.velocity[axis], 1e-12);
        close(
            receipt.exchange.hydrodynamic.solid_impulse[axis],
            0.0,
            1e-11,
        );
    }
}

#[test]
fn native_streaming_body_orbital_and_spin_momenta_close_over_moving_steps() {
    let mut grid = Grid::uniform(72, 72, 0.9);
    let initial = body([32.25, 33.125], [0.017, -0.008], 0.006);
    let mut disk = MovingDisk2::new(&grid, 2.4, initial, -0.2);
    let mut expected = totals(&grid, initial);
    for tick in 0..16 {
        let impulse = if tick < 8 {
            [0.002, -0.001]
        } else {
            [-0.001, 0.0015]
        };
        let couple = if tick % 2 == 0 { 0.001 } else { -0.0004 };
        let before = disk.body();
        let center = [
            before.center[0] + 0.5 * before.velocity[0],
            before.center[1] + 0.5 * before.velocity[1],
        ];
        expected[1] += impulse[0];
        expected[2] += impulse[1];
        let angular_impulse = center[0] * impulse[1] - center[1] * impulse[0] + couple;
        expected[3] += angular_impulse;
        let receipt = disk.step(&mut grid, impulse, couple);
        assert_eq!(receipt.collision_center, center);
        close(
            receipt.external_angular_impulse_about_origin,
            angular_impulse,
            1e-16,
        );
        let actual = totals(&grid, disk.body());
        close(actual[0], expected[0], 5e-9);
        for quantity in 1..4 {
            close(actual[quantity], expected[quantity], 2e-9);
        }
        for residual in receipt.exchange.momentum_residual {
            close(residual, 0.0, 2e-12);
        }
        // The domain is large enough that the localized perturbation has not
        // reached a periodic cut; no outer-boundary angular flux is omitted.
    }
    assert!(disk.body().center[0] > initial.center[0]);
    assert_ne!(disk.angle(), -0.2);
    assert_ne!(disk.body().angular_velocity, initial.angular_velocity);
}

#[test]
fn moving_geometry_and_native_fluid_replay_bit_for_bit() {
    let mut first_grid = Grid::uniform(28, 24, 0.85);
    let initial = body([12.25, 11.5], [0.022, 0.013], -0.008);
    let mut first = MovingDisk2::new(&first_grid, 2.2, initial, 0.6);
    let mut second_grid = first_grid.clone();
    let mut second = first.clone();
    for tick in 0..32 {
        let impulse = [0.001 * (tick % 3) as f64, -0.0007];
        let a = first.step(&mut first_grid, impulse, 0.0005);
        let b = second.step(&mut second_grid, impulse, 0.0005);
        assert_eq!(a, b);
        assert_eq!(first.body(), second.body());
        assert_eq!(first.angle(), second.angle());
        assert_eq!(first_grid.f, second_grid.f);
    }
}

#[test]
fn late_endpoint_refusal_preserves_body_fluid_and_scratch_recovery() {
    let mut grid = Grid::uniform(20, 14, 0.8);
    let initial = PartialSaturationBody2 {
        mass: 1000.0,
        inertia: 800.0,
        center: [16.97, 7.0],
        velocity: [0.0; 2],
        angular_velocity: 0.0,
    };
    let mut disk = MovingDisk2::new(&grid, 1.5, initial, 0.0);
    let mut control_grid = grid.clone();
    let mut control = disk.clone();
    let before_populations = grid.f.clone();
    // The midpoint is admitted and the native inertia solve succeeds below
    // the low-Mach cap, but the resulting final half-drift leaves the halo.
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = disk.step(&mut grid, [160.0, 0.0], 0.0);
        }))
        .is_err()
    );
    assert_eq!(disk.body(), initial);
    assert_eq!(disk.angle(), 0.0);
    assert_eq!(grid.f, before_populations);
    let actual = disk.step(&mut grid, [0.0; 2], 0.0);
    let expected = control.step(&mut control_grid, [0.0; 2], 0.0);
    assert_eq!(actual, expected);
    assert_eq!(disk.body(), control.body());
    assert_eq!(grid.f, control_grid.f);

    let wall = grid.idx(16, 7);
    grid.flags[wall] = Cell::Wall;
    let body_before = disk.body();
    let fluid_before = grid.f.clone();
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = disk.step(&mut grid, [0.0; 2], 0.0);
        }))
        .is_err()
    );
    assert_eq!(disk.body(), body_before);
    assert_eq!(grid.f, fluid_before);

    assert!(std::panic::catch_unwind(|| disk_cell_fraction(f64::INFINITY, [0.0; 2])).is_err());
    assert!(std::panic::catch_unwind(|| disk_cell_fraction(0.0, [0.0; 2])).is_err());
}
