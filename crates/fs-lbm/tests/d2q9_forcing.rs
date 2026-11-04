//! Physical momentum histories on the original dense collide/stream owner.

use fs_exec::CancelGate;
use fs_lbm::core2::{Cell, CollisionModel2, ForceRebindError2, ForceStepError2, Grid};
use fs_lbm::{CS2, Q, equilibrium};

const DIRECTIONS: [[f64; 2]; Q] = [
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

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-14 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

// Independent quadrature of density, PHYSICAL momentum, and the complete
// physical second moment. Equilibrium replacement cannot pass the stress test.
fn physical(grid: &Grid, cell: usize) -> [f64; 6] {
    let f = grid.f[cell];
    let rho = f.iter().sum::<f64>();
    let force: [f64; 2] = core::array::from_fn(|d| rho * grid.g[d] + grid.fext[cell][d]);
    let momentum: [f64; 2] = core::array::from_fn(|d| {
        f.iter().zip(DIRECTIONS).map(|(p, e)| p * e[d]).sum::<f64>() + 0.5 * force[d]
    });
    let u = momentum.map(|p| p / rho);
    let second = |a: usize, b: usize| {
        f.iter()
            .zip(DIRECTIONS)
            .map(|(p, e)| p * e[a] * e[b])
            .sum::<f64>()
            + 0.5 * (u[a] * force[b] + force[a] * u[b])
    };
    [
        rho,
        momentum[0],
        momentum[1],
        second(0, 0),
        second(0, 1),
        second(1, 1),
    ]
}

fn unchanged(actual: &Grid, before: &Grid) {
    assert_eq!(actual.flags, before.flags);
    assert_eq!(actual.tau, before.tau);
    assert_eq!(actual.g.map(f64::to_bits), before.g.map(f64::to_bits));
    for (a, b) in actual.f.iter().zip(&before.f) {
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    }
    for (a, b) in actual.fext.iter().zip(&before.fext) {
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    }
}

#[test]
fn rebind_preserves_nonuniform_density_momentum_and_full_physical_stress() {
    let mut grid = Grid::uniform(4, 3, 0.8);
    grid.flags[0] = Cell::Wall;
    grid.flags[1] = Cell::Gas;
    grid.flags[2] = Cell::Interface;
    for i in 0..grid.f.len() {
        grid.f[i] = equilibrium(0.6 + 0.11 * i as f64, 0.024, -0.019);
        // Nonzero xx stress plus an independent xy stress, both mass-free.
        grid.f[i][1] += 0.001;
        grid.f[i][3] += 0.001;
        grid.f[i][0] -= 0.002;
        grid.f[i][5] += 0.0003;
        grid.f[i][7] += 0.0003;
        grid.f[i][6] -= 0.0003;
        grid.f[i][8] -= 0.0003;
        grid.tau[i] += 0.01 * i as f64;
    }
    let initial = grid.clone();
    let observations: Vec<_> = (0..grid.f.len()).map(|i| physical(&grid, i)).collect();
    let mut scratch = Vec::new();
    for scale in [1.0, -1.7, 0.0] {
        let external: Vec<_> = (0..grid.f.len())
            .map(|i| [scale * 0.0004 * i as f64, scale * -0.0007])
            .collect();
        grid.rebind_forcing([scale * 0.003, scale * -0.002], &external, &mut scratch)
            .unwrap();
        for (i, expected) in observations.iter().enumerate().skip(2) {
            for (actual, expected) in physical(&grid, i).into_iter().zip(expected) {
                near(actual, *expected);
            }
        }
        assert_eq!(grid.f[..2], initial.f[..2], "inactive populations changed");
        assert_eq!(grid.flags, initial.flags);
        assert_eq!(grid.tau, initial.tau);
        assert_eq!(grid.fext, external);

        let accepted = grid.clone();
        let scratch_before = scratch.clone();
        grid.rebind_forcing(grid.g, &external, &mut scratch)
            .unwrap();
        unchanged(&grid, &accepted);
        assert_eq!(scratch, scratch_before, "identical force must be a no-op");
    }

    // Changing only gravity retains the independently supplied local load.
    let external = vec![[0.002, -0.001]; grid.f.len()];
    grid.rebind_forcing([0.0; 2], &external, &mut scratch)
        .unwrap();
    grid.rebind_gravity([-0.002, 0.004], &mut scratch).unwrap();
    assert_eq!(grid.fext, external);
    for (i, expected) in observations.iter().enumerate().skip(2) {
        for (actual, expected) in physical(&grid, i).into_iter().zip(expected) {
            near(actual, *expected);
        }
    }
}

#[test]
fn switching_and_reversing_forces_delivers_only_the_elapsed_impulse() {
    for collision in [
        CollisionModel2::Bgk,
        CollisionModel2::Regularized,
        CollisionModel2::RecursiveRegularized,
        CollisionModel2::Trt { magic: 3.0 / 16.0 },
        CollisionModel2::CentralMoment,
    ] {
        for rho in [0.55, 1.9] {
            let mut grid = Grid::uniform(3, 2, 0.83);
            grid.collision = collision;
            let mut expected_velocity = [0.023, -0.017];
            grid.f
                .fill(equilibrium(rho, expected_velocity[0], expected_velocity[1]));
            let mut scratch = Vec::new();
            for (g, local) in [
                ([0.002, -0.001], [0.0007, 0.0002]),
                ([-0.003, 0.002], [-0.0004, 0.0003]),
                ([0.0, 0.001], [0.0005, -0.0001]),
                ([0.0; 2], [0.0; 2]),
            ] {
                grid.rebind_forcing(g, &vec![local; grid.f.len()], &mut scratch)
                    .unwrap();
                for i in 0..grid.f.len() {
                    let before_step = grid.moments(i);
                    near(before_step.rho, rho);
                    for d in 0..2 {
                        near(before_step.u[d], expected_velocity[d]);
                    }
                }
                grid.step(&mut scratch);
                for d in 0..2 {
                    expected_velocity[d] += g[d] + local[d] / rho;
                }
                for i in 0..grid.f.len() {
                    let after_step = grid.moments(i);
                    near(after_step.rho, rho);
                    for d in 0..2 {
                        near(after_step.u[d], expected_velocity[d]);
                    }
                }
            }
        }
    }
}

#[test]
fn nonuniform_periodic_flow_closes_total_momentum_against_actual_body_force() {
    let mut grid = Grid::uniform(5, 4, 0.9);
    for (i, f) in grid.f.iter_mut().enumerate() {
        *f = equilibrium(0.8 + 0.02 * i as f64, 0.002 * (i % 3) as f64, -0.01);
    }
    let external: Vec<_> = (0..grid.f.len())
        .map(|i| [0.0001 * (i % 4) as f64, -0.0002 * (i % 3) as f64])
        .collect();
    let mut scratch = Vec::new();
    grid.rebind_forcing([0.001, -0.0007], &external, &mut scratch)
        .unwrap();
    let mut before = [0.0; 3];
    let mut supplied = [0.0; 2];
    for i in 0..grid.f.len() {
        let state = physical(&grid, i);
        for d in 0..3 {
            before[d] += state[d];
        }
        for d in 0..2 {
            supplied[d] += state[0] * grid.g[d] + external[i][d];
        }
    }
    grid.step(&mut scratch);
    let mut after = [0.0; 3];
    for i in 0..grid.f.len() {
        let state = physical(&grid, i);
        for d in 0..3 {
            after[d] += state[d];
        }
    }
    near(after[0], before[0]);
    for d in 0..2 {
        near(after[d + 1] - before[d + 1], supplied[d]);
    }
}

#[test]
fn endpoint_force_ramps_and_reversals_integrate_the_physical_impulse() {
    let gate = CancelGate::new();
    for collision in [
        CollisionModel2::Bgk,
        CollisionModel2::Regularized,
        CollisionModel2::RecursiveRegularized,
        CollisionModel2::Trt { magic: 3.0 / 16.0 },
        CollisionModel2::CentralMoment,
    ] {
        for (rho, tau) in [(0.55, 0.61), (1.9, 1.37)] {
            let mut grid = Grid::uniform(3, 2, tau);
            grid.collision = collision;
            let mut expected_velocity = [0.023, -0.017];
            grid.f
                .fill(equilibrium(rho, expected_velocity[0], expected_velocity[1]));
            let mut post = Vec::new();
            let mut next = Vec::new();
            for (gravity, local) in [
                ([0.002, -0.001], [0.0007, 0.0002]),
                ([0.004, -0.002], [0.0014, 0.0004]),
                ([-0.003, 0.002], [-0.0004, 0.0003]),
                ([0.0; 2], [0.0; 2]),
            ] {
                for d in 0..2 {
                    expected_velocity[d] +=
                        0.5 * (grid.g[d] + gravity[d] + (grid.fext[0][d] + local[d]) / rho);
                }
                grid.step_with_endpoint_forcing(
                    gravity,
                    &vec![local; grid.f.len()],
                    &mut post,
                    &mut next,
                    &gate,
                )
                .unwrap();
                for cell in 0..grid.f.len() {
                    let state = physical(&grid, cell);
                    near(state[0], rho);
                    for d in 0..2 {
                        near(state[d + 1] / state[0], expected_velocity[d]);
                    }
                }
            }
        }
    }
}

#[test]
fn endpoint_forces_close_spatial_momentum_with_native_wall_reaction() {
    // Cross a cancellation chunk boundary; retain nonuniform density and a
    // stationary internal wall so neither uniform-state nor periodic-only
    // momentum accounting can satisfy this balance accidentally.
    let mut grid = Grid::uniform(18, 17, 0.9);
    grid.flags[..18].fill(Cell::Wall);
    for (i, f) in grid.f.iter_mut().enumerate() {
        *f = equilibrium(0.8 + 0.01 * (i % 23) as f64, 0.003 * (i % 3) as f64, -0.01);
    }
    let initial_external: Vec<_> = (0..grid.f.len())
        .map(|i| [0.0001 * (i % 4) as f64, -0.0002 * (i % 3) as f64])
        .collect();
    let mut post = Vec::new();
    let mut next = Vec::new();
    grid.rebind_forcing([0.001, -0.0007], &initial_external, &mut post)
        .unwrap();
    let before = grid.clone();
    let mut native = grid.clone();
    native.collide_into(&mut post);
    let wall_mask: Vec<_> = grid.flags.iter().map(|flag| *flag == Cell::Wall).collect();
    let reaction = native.stream_from_with_wall_momentum(&post, &wall_mask);
    assert!(reaction.measured_links > 0);

    let endpoint_gravity = [-0.0004, 0.0013];
    let endpoint_external: Vec<_> = (0..grid.f.len())
        .map(|i| [-0.0003 * (i % 5) as f64, 0.0001 * (i % 7) as f64])
        .collect();
    grid.step_with_endpoint_forcing(
        endpoint_gravity,
        &endpoint_external,
        &mut post,
        &mut next,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        grid.f, native.f,
        "endpoint source must retain the native pull"
    );
    assert_eq!(grid.f[..18], before.f[..18]);

    let mut mass_change = 0.0;
    let mut momentum_change = [0.0; 2];
    let mut supplied = [0.0; 2];
    let mut density_transport = 0.0;
    for cell in 18..grid.f.len() {
        let previous = physical(&before, cell);
        let current = physical(&grid, cell);
        mass_change += current[0] - previous[0];
        density_transport += (current[0] - previous[0]).abs();
        for d in 0..2 {
            momentum_change[d] += current[d + 1] - previous[d + 1];
            supplied[d] += 0.5
                * (previous[0] * before.g[d]
                    + initial_external[cell][d]
                    + current[0] * endpoint_gravity[d]
                    + endpoint_external[cell][d]);
        }
    }
    assert!(density_transport > 0.01, "fixture must transport density");
    near(mass_change, 0.0);
    for d in 0..2 {
        near(momentum_change[d] + reaction.wall_impulse[d], supplied[d]);
    }
}

#[test]
fn unchanged_endpoint_forces_preserve_every_native_collision_bit() {
    let gate = CancelGate::new();
    for collision in [
        CollisionModel2::Bgk,
        CollisionModel2::Regularized,
        CollisionModel2::RecursiveRegularized,
        CollisionModel2::Trt { magic: 3.0 / 16.0 },
        CollisionModel2::CentralMoment,
    ] {
        let mut grid = Grid::uniform(19, 15, 0.8);
        grid.collision = collision;
        grid.periodic_x = false;
        grid.flags[20] = Cell::Wall;
        grid.flags[21] = Cell::Gas;
        grid.flags[22] = Cell::Interface;
        for (i, f) in grid.f.iter_mut().enumerate() {
            *f = equilibrium(0.8 + 0.02 * (i % 11) as f64, 0.01, -0.012);
            grid.tau[i] += 0.02 * (i % 7) as f64;
        }
        let external = vec![[0.0002, -0.0001]; grid.f.len()];
        let mut post = Vec::new();
        let mut next = Vec::new();
        grid.rebind_forcing([0.001, -0.002], &external, &mut post)
            .unwrap();
        let mut native = grid.clone();
        for _ in 0..3 {
            native.step(&mut post);
            grid.step_with_endpoint_forcing(grid.g, &external, &mut post, &mut next, &gate)
                .unwrap();
            unchanged(&grid, &native);
        }
    }
}

#[test]
fn refused_endpoint_steps_publish_neither_populations_nor_forces() {
    let initial = Grid::uniform(3, 2, 0.8);
    let mut post = Vec::new();
    let mut next = Vec::new();
    for case in 0..10 {
        let mut grid = initial.clone();
        let mut gravity = [0.001, -0.002];
        let mut external = vec![[0.0002, 0.0001]; grid.f.len()];
        let gate = CancelGate::new();
        match case {
            0 => {
                external.pop();
            }
            1 => gravity[1] = f64::NAN,
            2 => external[5][0] = f64::NAN,
            3 => grid.f[5] = [0.0; Q],
            4 => grid.tau[5] = 0.5,
            5 => {
                grid.tau.pop();
            }
            6 => grid.collision = CollisionModel2::Trt { magic: 0.0 },
            7 => {
                grid.f.fill(equilibrium(2.0, 0.0, 0.0));
                gravity[0] = f64::MAX;
            }
            // The inputs and initial physical velocity are finite; collision
            // overflows. This must refuse after proposal work, before publish.
            8 => grid.fext[5][0] = 1e200,
            _ => gate.request(),
        }
        let before = grid.clone();
        assert!(matches!(
            grid.step_with_endpoint_forcing(gravity, &external, &mut post, &mut next, &gate),
            Err(ForceStepError2::Invalid(_)
                | ForceStepError2::Unphysical { .. }
                | ForceStepError2::Cancelled)
        ));
        unchanged(&grid, &before);
    }
}

#[test]
fn rejected_force_changes_do_not_publish_early_cells() {
    let initial = Grid::uniform(3, 2, 0.8);
    let mut scratch = Vec::new();
    for case in 0..6 {
        let mut grid = initial.clone();
        let mut g = [0.001, -0.002];
        let mut external = vec![[0.0002, 0.0001]; grid.f.len()];
        match case {
            0 => {
                external.pop();
            }
            1 => g[1] = f64::NAN,
            2 => external[5][0] = f64::NAN,
            3 => grid.f[5] = [0.0; Q],
            4 => {
                grid.f[5] = equilibrium(2.0, 0.0, 0.0);
                g[0] = f64::MAX;
            }
            _ => grid.fext[5][1] = f64::INFINITY,
        }
        let before = grid.clone();
        assert!(matches!(
            grid.rebind_forcing(g, &external, &mut scratch),
            Err(ForceRebindError2::Invalid(_) | ForceRebindError2::Unphysical { .. })
        ));
        unchanged(&grid, &before);
    }
}

#[test]
fn ordinary_dense_initializers_start_at_physical_rest() {
    let g = 0.002;
    let channel = fs_lbm::Lbm::channel(3, 6, 0.8, g);
    for y in 0..6 {
        let (ux, uy) = channel.velocity(1, y);
        near(ux, 0.0);
        near(uy, 0.0);
    }
    let dam = fs_lbm::freesurface::dam_break(12, 12, 3, g, 0.0, fs_lbm::ContactModel::Neutral);
    let refined = fs_lbm::RefinedChannel::new(3, 4, 8, 0.8, g);
    let (rheology, steps, _) =
        fs_lbm::rheology::channel_flow(3, 6, g, fs_lbm::Rheology::Newtonian { nu: 0.1 }, 0);
    assert_eq!(steps, 0);
    for grid in [&dam.grid, &refined.coarse, &refined.fine, &rheology] {
        for i in 0..grid.f.len() {
            if matches!(grid.flags[i], Cell::Fluid | Cell::Interface) {
                let state = physical(grid, i);
                near(state[1], 0.0);
                near(state[2], 0.0);
                near(state[3], state[0] * CS2);
                near(state[4], 0.0);
                near(state[5], state[0] * CS2);
            }
        }
    }
}
