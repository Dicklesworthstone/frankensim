//! G0/G1/G3/G4 checks of reciprocal immersed-body momentum coupling.

use fs_lbm::core2::{CollisionModel2, PartialSaturationBody2};
use fs_lbm::{Cell, Grid, Q, equilibrium};

const EX: [f64; Q] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const EY: [f64; Q] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

fn fluid_momentum(populations: &[[f64; Q]], nx: usize, origin: [f64; 2]) -> [f64; 3] {
    let mut total = [0.0; 3];
    for (cell, f) in populations.iter().enumerate() {
        let x = (cell % nx) as f64 - origin[0];
        let y = (cell / nx) as f64 - origin[1];
        for q in 0..Q {
            total[0] += EX[q] * f[q];
            total[1] += EY[q] * f[q];
            total[2] += (x * EY[q] - y * EX[q]) * f[q];
        }
    }
    total
}

fn total_momentum(
    populations: &[[f64; Q]],
    nx: usize,
    body: PartialSaturationBody2,
    origin: [f64; 2],
) -> [f64; 3] {
    let mut p = fluid_momentum(populations, nx, origin);
    p[0] += body.mass * body.velocity[0];
    p[1] += body.mass * body.velocity[1];
    p[2] += body.inertia * body.angular_velocity
        + body.mass
            * ((body.center[0] - origin[0]) * body.velocity[1]
                - (body.center[1] - origin[1]) * body.velocity[0]);
    p
}

fn kinetic_energy(populations: &[[f64; Q]], body: PartialSaturationBody2) -> f64 {
    let mut energy = 0.5
        * (body.mass * (body.velocity[0].powi(2) + body.velocity[1].powi(2))
            + body.inertia * body.angular_velocity.powi(2));
    for f in populations {
        let rho: f64 = f.iter().sum();
        let px: f64 = (0..Q).map(|q| EX[q] * f[q]).sum();
        let py: f64 = (0..Q).map(|q| EY[q] * f[q]).sum();
        energy += 0.5 * (px * px + py * py) / rho;
    }
    energy
}

fn body() -> PartialSaturationBody2 {
    PartialSaturationBody2 {
        mass: 0.7,
        inertia: 0.9,
        center: [1.7, 1.8],
        velocity: [-0.01, 0.015],
        angular_velocity: 0.012,
    }
}

#[test]
fn g1_off_center_cell_uses_the_coupled_point_mobility() {
    let mut grid = Grid::uniform(6, 5, 0.9);
    let cell = grid.idx(3, 1);
    let rho = 1.1;
    let fluid_velocity = [0.04, -0.02];
    grid.f[cell] = equilibrium(rho, fluid_velocity[0], fluid_velocity[1]);
    let mut fractions = vec![0.0; grid.f.len()];
    fractions[cell] = 0.7;
    let before = body();
    let mut moving = before;

    // Independently eliminate body velocities in favor of the TWO unknown
    // components of impulse J at this one cell. Point mobility is
    // 1/m + (r cross)^T (1/I) (r cross), including its off-diagonal terms.
    let r = [3.0 - before.center[0], 1.0 - before.center[1]];
    let gamma = rho * 0.7 * (0.9 - 0.5) / (1.0 - 0.7 + 0.9 - 0.5);
    let slip = [
        fluid_velocity[0] - before.velocity[0] + before.angular_velocity * r[1],
        fluid_velocity[1] - before.velocity[1] - before.angular_velocity * r[0],
    ];
    let a = 1.0 + 0.5 * gamma * (1.0 / before.mass + r[1] * r[1] / before.inertia);
    let b = -0.5 * gamma * r[0] * r[1] / before.inertia;
    let d = 1.0 + 0.5 * gamma * (1.0 / before.mass + r[0] * r[0] / before.inertia);
    let determinant = a * d - b * b;
    let expected_impulse = [
        gamma * (d * slip[0] - b * slip[1]) / determinant,
        gamma * (a * slip[1] - b * slip[0]) / determinant,
    ];
    let expected_torque = r[0] * expected_impulse[1] - r[1] * expected_impulse[0];
    let mut post = Vec::new();
    let receipt = grid.collide_into_with_partial_saturation_body(
        &mut post,
        &fractions,
        &mut moving,
        [0.0; 2],
        0.0,
    );
    for (axis, &impulse) in expected_impulse.iter().enumerate() {
        assert!((receipt.hydrodynamic.solid_impulse[axis] - impulse).abs() < 2e-15);
        assert!(
            (moving.velocity[axis] - before.velocity[axis] - impulse / before.mass).abs() < 2e-15
        );
    }
    assert!((receipt.hydrodynamic.solid_angular_impulse - expected_torque).abs() < 2e-15);
    assert!(
        (moving.angular_velocity - before.angular_velocity - expected_torque / before.inertia)
            .abs()
            < 2e-15
    );
    assert!((receipt.kinetic_energy_change - receipt.hydrodynamic.solid_work).abs() < 2e-15);
}

#[test]
fn g0_actual_streamed_populations_close_linear_orbital_and_spin_momentum() {
    let mut grid = Grid::uniform(9, 9, 0.8);
    let mut fractions = vec![0.0; grid.f.len()];
    for (x, y, fraction, rho, velocity, tau) in [
        (3, 3, 0.25, 1.03, [0.03, -0.01], 0.6),
        (4, 5, 0.8, 0.98, [-0.02, 0.04], 1.3),
        (5, 4, 1.0, 1.1, [0.01, 0.025], 0.9),
    ] {
        let cell = grid.idx(x, y);
        fractions[cell] = fraction;
        grid.tau[cell] = tau;
        grid.f[cell] = equilibrium(rho, velocity[0], velocity[1]);
    }
    let mut moving = PartialSaturationBody2 {
        center: [3.7, 4.2],
        ..body()
    };
    let origin = [-0.2, 0.4];
    let before = total_momentum(&grid.f, grid.nx, moving, origin);
    let mass_before: f64 = grid.f.iter().flatten().sum();
    let external = [0.003, -0.002];
    let external_couple = 0.001;
    let mut scratch = Vec::new();
    let receipt = grid.step_with_partial_saturation_body(
        &mut scratch,
        &fractions,
        &mut moving,
        external,
        external_couple,
    );
    let after = total_momentum(&grid.f, grid.nx, moving, origin);
    let angular_external = external_couple + (moving.center[0] - origin[0]) * external[1]
        - (moving.center[1] - origin[1]) * external[0];
    // The perturbed patch is interior. Uniform resting populations at the
    // periodic seams contribute equal opposite angular transport this step.
    for (axis, load) in [external[0], external[1], angular_external]
        .into_iter()
        .enumerate()
    {
        assert!((after[axis] - before[axis] - load).abs() < 2e-13);
        assert!(receipt.momentum_residual[axis].abs() < 2e-15);
    }
    let mass_after: f64 = grid.f.iter().flatten().sum();
    assert!((mass_after - mass_before).abs() < 2e-12);
    assert!(
        (receipt.kinetic_energy_change - receipt.hydrodynamic.solid_work - receipt.external_work)
            .abs()
            < 2e-15
    );
}

#[test]
fn g0_light_translating_and_spinning_body_is_passive_and_replays() {
    let mut grid = Grid::uniform(5, 5, 0.8);
    let mut fractions = vec![0.0; grid.f.len()];
    let mut moving = PartialSaturationBody2 {
        mass: 1e-6,
        inertia: 2e-6,
        center: [2.0, 2.0],
        velocity: [-0.015, 0.02],
        angular_velocity: -0.025,
    };
    for (x, y) in [(1, 2), (3, 2), (2, 1), (2, 3)] {
        let cell = grid.idx(x, y);
        fractions[cell] = 1.0;
        let u = [
            0.02 - 0.015 * (y as f64 - 2.0),
            -0.01 + 0.015 * (x as f64 - 2.0),
        ];
        grid.f[cell] = equilibrium(1.0, u[0], u[1]);
    }
    let mut replay = grid.clone();
    let mut replay_body = moving;
    let p0 = total_momentum(&grid.f, grid.nx, moving, [-1.0, 0.5]);
    let mut prior_energy = kinetic_energy(&grid.f, moving);
    let mut post = Vec::new();
    let mut replay_post = Vec::new();
    for _ in 0..128 {
        let receipt = grid.collide_into_with_partial_saturation_body(
            &mut post,
            &fractions,
            &mut moving,
            [0.0; 2],
            0.0,
        );
        let other = replay.collide_into_with_partial_saturation_body(
            &mut replay_post,
            &fractions,
            &mut replay_body,
            [0.0; 2],
            0.0,
        );
        assert_eq!(receipt, other);
        assert_eq!(moving, replay_body);
        assert_eq!(post, replay_post);
        let energy = kinetic_energy(&post, moving);
        assert!(
            energy <= prior_energy + 2e-16,
            "coupling created energy: {prior_energy} -> {energy}"
        );
        prior_energy = energy;
        // Repeated collisions isolate the passive fluid/body exchange from
        // compressible transport and its pressure-work/energy accounting.
        std::mem::swap(&mut grid.f, &mut post);
        std::mem::swap(&mut replay.f, &mut replay_post);
    }
    let p1 = total_momentum(&grid.f, grid.nx, moving, [-1.0, 0.5]);
    for axis in 0..3 {
        assert!((p1[axis] - p0[axis]).abs() < 2e-12);
    }
}

#[test]
fn g3_common_rigid_motion_has_no_spurious_hydrodynamic_kick() {
    let mut grid = Grid::uniform(6, 5, 0.8);
    let mut fractions = vec![0.0; grid.f.len()];
    let mut moving = body();
    let before = moving;
    for (x, y, fraction) in [(1, 1, 0.25), (3, 1, 0.8), (3, 3, 1.0)] {
        let cell = grid.idx(x, y);
        fractions[cell] = fraction;
        let u = [
            moving.velocity[0] - moving.angular_velocity * (y as f64 - moving.center[1]),
            moving.velocity[1] + moving.angular_velocity * (x as f64 - moving.center[0]),
        ];
        grid.f[cell] = equilibrium(1.1, u[0], u[1]);
    }
    let receipt = grid.collide_into_with_partial_saturation_body(
        &mut Vec::new(),
        &fractions,
        &mut moving,
        [0.0; 2],
        0.0,
    );
    for axis in 0..2 {
        assert!((moving.velocity[axis] - before.velocity[axis]).abs() < 2e-15);
        assert!(receipt.hydrodynamic.solid_impulse[axis].abs() < 2e-15);
    }
    assert!((moving.angular_velocity - before.angular_velocity).abs() < 2e-15);
    assert!(receipt.hydrodynamic.solid_angular_impulse.abs() < 2e-15);
}

#[test]
fn g0_zero_fraction_keeps_ordinary_fluid_and_free_body_impulse() {
    let grid = Grid::uniform(3, 3, 0.8);
    let mut ordinary = Vec::new();
    grid.collide_into(&mut ordinary);
    let before = body();
    let mut moving = before;
    let external = [0.002, -0.003];
    let couple = 0.004;
    let mut post = Vec::new();
    let receipt = grid.collide_into_with_partial_saturation_body(
        &mut post,
        &[0.0; 9],
        &mut moving,
        external,
        couple,
    );
    assert_eq!(ordinary, post);
    assert_eq!(receipt.hydrodynamic.coupled_cells, 0);
    for axis in 0..2 {
        assert!(
            (moving.velocity[axis] - before.velocity[axis] - external[axis] / before.mass).abs()
                < 2e-15
        );
    }
    assert!(
        (moving.angular_velocity - before.angular_velocity - couple / before.inertia).abs() < 2e-15
    );
    assert!((receipt.kinetic_energy_change - receipt.external_work).abs() < 2e-15);
}

fn body_bits(body: PartialSaturationBody2) -> [u64; 7] {
    [
        body.mass,
        body.inertia,
        body.center[0],
        body.center[1],
        body.velocity[0],
        body.velocity[1],
        body.angular_velocity,
    ]
    .map(f64::to_bits)
}

#[test]
fn g4_invalid_body_geometry_and_late_velocity_refusal_publish_nothing() {
    for case in 0..12 {
        let mut grid = Grid::uniform(5, 5, 0.8);
        let mut moving = body();
        let mut fractions = vec![0.0; grid.f.len()];
        fractions[grid.idx(2, 2)] = 0.7;
        let mut external = [0.0; 2];
        match case {
            0 => moving.mass = 0.0,
            1 => moving.inertia = -1.0,
            2 => moving.mass = f64::NAN,
            3 => moving.center[0] = f64::INFINITY,
            4 => fractions[12] = -0.1,
            5 => fractions[12] = 1.1,
            6 => fractions[12] = f64::NAN,
            7 => grid.g = [0.001, 0.0],
            8 => grid.collision = CollisionModel2::Regularized,
            9 => grid.flags[0] = Cell::Interface,
            10 => external = [1e4, 0.0],
            _ => moving.angular_velocity = 1e4,
        }
        let before_populations = grid.f.clone();
        let before_body = body_bits(moving);
        let mut scratch = vec![[123.0; Q]; 2];
        let before_scratch = scratch.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = grid.step_with_partial_saturation_body(
                &mut scratch,
                &fractions,
                &mut moving,
                external,
                0.0,
            );
        }));
        assert!(result.is_err(), "case {case} was admitted");
        assert_eq!(grid.f, before_populations);
        assert_eq!(body_bits(moving), before_body);
        assert_eq!(scratch, before_scratch);
    }
}
