use fs_lbm::core2::frame::{FrameForcingError2, RigidFrame2};
use fs_lbm::{Grid, equilibrium};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn changing_frame_force_preserves_current_physical_momentum_and_stress() {
    let mut grid = Grid::uniform(1, 1, 0.8);
    grid.f[0] = equilibrium(1.7, 0.03, -0.02);
    grid.g = [0.001, -0.002];
    grid.fext[0] = [0.002, 0.001];
    let before = grid.moments(0);
    let physical_second = |grid: &Grid| {
        let directions = [
            [0., 0.],
            [1., 0.],
            [0., 1.],
            [-1., 0.],
            [0., -1.],
            [1., 1.],
            [-1., 1.],
            [-1., -1.],
            [1., -1.],
        ];
        let m = grid.moments(0);
        let force: [f64; 2] = core::array::from_fn(|d| grid.g[d].mul_add(m.rho, grid.fext[0][d]));
        core::array::from_fn::<_, 2, _>(|a| {
            core::array::from_fn::<_, 2, _>(|b| {
                grid.f[0]
                    .iter()
                    .zip(directions)
                    .map(|(f, e)| f * e[a] * e[b])
                    .sum::<f64>()
                    + 0.5 * (m.u[a] * force[b] + m.u[b] * force[a])
            })
        })
    };
    let stress_before = physical_second(&grid);
    RigidFrame2 {
        origin: [2., -3.],
        angular_velocity: 0.004,
        angular_acceleration: -0.0001,
        origin_acceleration: [0.0002, -0.0003],
    }
    .prepare_step(&mut grid, [0.0, -0.0006])
    .unwrap();
    let after = grid.moments(0);
    near(after.rho, before.rho, 1e-15);
    let stress_after = physical_second(&grid);
    for a in 0..2 {
        near(after.u[a], before.u[a], 1e-16);
        for b in 0..2 {
            near(stress_after[a][b], stress_before[a][b], 2e-16);
        }
    }
}

#[test]
fn coriolis_source_rotates_momentum_without_adding_kinetic_energy() {
    let omega = 0.13;
    let frame = RigidFrame2 {
        angular_velocity: omega,
        ..RigidFrame2::default()
    };
    let mut grid = Grid::uniform(1, 1, 0.8);
    let initial = [0.03, -0.02];
    grid.f[0] = equilibrium(1.4, initial[0], initial[1]);
    let mut scratch = Vec::new();
    let mut expected = initial;
    for _ in 0..2000 {
        frame.prepare_step(&mut grid, [0.0; 2]).unwrap();
        grid.step(&mut scratch);
        expected = [
            ((1.0 - omega * omega) * expected[0] + 2.0 * omega * expected[1])
                / (1.0 + omega * omega),
            ((1.0 - omega * omega) * expected[1] - 2.0 * omega * expected[0])
                / (1.0 + omega * omega),
        ];
    }
    let actual = grid.moments(0);
    near(actual.rho, 1.4, 2e-12);
    near(actual.u[0], expected[0], 2e-13);
    near(actual.u[1], expected[1], 2e-13);
    near(
        actual.u[0] * actual.u[0] + actual.u[1] * actual.u[1],
        initial[0] * initial[0] + initial[1] * initial[1],
        2e-14,
    );
}

#[test]
fn angular_acceleration_applies_the_opposing_euler_impulse() {
    let mut grid = Grid::uniform(1, 1, 0.9);
    let frame = RigidFrame2 {
        origin: [-2., 3.],
        angular_acceleration: 0.001,
        ..RigidFrame2::default()
    };
    frame.prepare_step(&mut grid, [0.0; 2]).unwrap();
    let start = grid.moments(0).u;
    near(start[0], 0.0, 1e-16);
    near(start[1], 0.0, 1e-16);
    grid.step(&mut Vec::new());
    let end = grid.moments(0).u;
    near(end[0], -0.003, 1e-16);
    near(end[1], -0.002, 1e-16);
    near(2.0 * end[1] + 3.0 * end[0], -0.013, 1e-15);
}

#[test]
fn centrifugal_and_coriolis_forces_share_the_same_midpoint_velocity() {
    let mut grid = Grid::uniform(1, 1, 0.75);
    let omega = 0.02;
    let frame = RigidFrame2 {
        origin: [-2.0, 0.0],
        angular_velocity: omega,
        ..RigidFrame2::default()
    };
    frame.prepare_step(&mut grid, [0.0; 2]).unwrap();
    grid.step(&mut Vec::new());
    let velocity = grid.moments(0).u;
    near(
        velocity[0],
        2.0 * omega * omega / (1.0 + omega * omega),
        1e-16,
    );
    near(
        velocity[1],
        -2.0 * omega * omega * omega / (1.0 + omega * omega),
        1e-16,
    );
}

#[test]
fn freely_falling_frame_cancels_uniform_gravity() {
    let mut grid = Grid::uniform(3, 2, 0.8);
    let before = grid.f.clone();
    let gravity = [0.003, -0.007];
    RigidFrame2 {
        origin_acceleration: gravity,
        ..RigidFrame2::default()
    }
    .prepare_step(&mut grid, gravity)
    .unwrap();
    assert_eq!(grid.f, before);
    assert!(grid.fext.iter().all(|force| *force == [0.0; 2]));
    grid.step(&mut Vec::new());
    for cell in 0..6 {
        assert_eq!(grid.moments(cell).u, [0.0; 2]);
    }
}

#[test]
fn invalid_late_cell_does_not_partially_reframe_the_grid() {
    let mut grid = Grid::uniform(3, 2, 0.8);
    grid.f[5] = [-1.0; 9];
    let before = grid.clone();
    assert_eq!(
        RigidFrame2::default().prepare_step(&mut grid, [0.001, 0.0]),
        Err(FrameForcingError2::InvalidCell { cell: 5 })
    );
    assert_eq!(grid.f, before.f);
    assert_eq!(grid.fext, before.fext);
    assert_eq!(grid.g, before.g);
    assert_eq!(
        RigidFrame2 {
            angular_velocity: f64::INFINITY,
            ..RigidFrame2::default()
        }
        .prepare_step(&mut grid, [0.0; 2]),
        Err(FrameForcingError2::InvalidMotion)
    );
}
