//! Native two-pressure D2Q9 flow: source-frame reconstruction, hydrostatics,
//! and pressure-head momentum transfer to an obstructing wall.

use fs_lbm::core2::{CollisionModel2, PressureProfileX2, PressureX2};
use fs_lbm::{CS2, Cell, Grid, Q, equilibrium};

const EX: [f64; Q] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const EY: [f64; Q] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

fn physical_stress(grid: &Grid, index: usize) -> [[f64; 2]; 2] {
    let moments = grid.moments(index);
    let eq = equilibrium(moments.rho, moments.u[0], moments.u[1]);
    let force: [f64; 2] =
        core::array::from_fn(|d| grid.g[d].mul_add(moments.rho, grid.fext[index][d]));
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            let velocity = [EX, EY];
            let raw = (0..Q)
                .map(|q| velocity[a][q] * velocity[b][q] * (grid.f[index][q] - eq[q]))
                .sum::<f64>();
            raw + 0.5 * (moments.u[a] * force[b] + moments.u[b] * force[a])
        })
    })
}

#[test]
fn prescribed_pressures_keep_the_local_guo_frame_and_original_wall_exchange() {
    let boundary = PressureX2::new(1.012, 0.997);
    assert_eq!(boundary.min_density(), 1.012);
    assert_eq!(boundary.max_density(), 0.997);
    let mut grid = Grid::uniform(9, 7, 0.83);
    grid.periodic_x = false;
    let mut forcing = vec![[0.0; 2]; grid.f.len()];
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let cell = grid.idx(x, y);
            grid.f[cell] = equilibrium(1.002 + 0.0001 * y as f64, 0.018, -0.003);
            // A nonzero, density- and momentum-conserving xx/shear stress.
            grid.f[cell][0] -= 0.0004;
            grid.f[cell][1] += 0.0002;
            grid.f[cell][3] += 0.0002;
            grid.f[cell][5] += 0.00002;
            grid.f[cell][6] -= 0.00002;
            grid.f[cell][7] += 0.00002;
            grid.f[cell][8] -= 0.00002;
            forcing[cell] = [1e-5 * (x + 1) as f64, -7e-6 * (x + y + 1) as f64];
        }
    }
    let wall = grid.idx(4, 3);
    grid.flags[wall] = Cell::Wall;
    grid.rebind_forcing([2e-5, -3e-5], &forcing, &mut Vec::new())
        .unwrap();
    let mut mask = vec![false; grid.f.len()];
    mask[wall] = true;
    let (mut scratch, mut measured_scratch) = (Vec::new(), Vec::new());
    let mut nonzero_stress = 0.0_f64;
    for _ in 0..4 {
        let mut measured = grid.clone();
        let mut original_stream = grid.clone();
        let expected = original_stream.step_with_wall_momentum(&mut Vec::new(), &mask);
        grid.step_pressure_x(&mut scratch, boundary);
        let observed =
            measured.step_pressure_x_with_wall_momentum(&mut measured_scratch, boundary, &mask);
        assert!(observed.measured_links > 0);
        assert_eq!(observed.measured_links, expected.measured_links);
        assert_eq!(
            observed.wall_impulse.map(f64::to_bits),
            expected.wall_impulse.map(f64::to_bits)
        );
        assert_eq!(scratch.len(), grid.f.len());
        for (plain, with_receipt) in grid.f.iter().zip(&measured.f) {
            assert_eq!(plain.map(f64::to_bits), with_receipt.map(f64::to_bits));
        }
        for y in 0..grid.ny {
            for (face_x, interior_x, rho) in [
                (0, 1, boundary.min_density()),
                (grid.nx - 1, grid.nx - 2, boundary.max_density()),
            ] {
                let face = grid.idx(face_x, y);
                let interior = grid.idx(interior_x, y);
                let actual = grid.moments(face);
                let donor = grid.moments(interior);
                assert!((actual.rho - rho).abs() < 2e-12);
                for (actual, donor) in actual.u.iter().zip(donor.u) {
                    assert!((actual - donor).abs() < 2e-12);
                }
                for (actual_row, donor_row) in physical_stress(&grid, face)
                    .iter()
                    .zip(physical_stress(&grid, interior))
                {
                    for (actual, donor) in actual_row.iter().zip(donor_row) {
                        nonzero_stress = nonzero_stress.max(donor.abs());
                        assert!((actual - donor).abs() < 2e-12);
                    }
                }
            }
            for x in 1..grid.nx - 1 {
                let index = grid.idx(x, y);
                assert_eq!(
                    grid.f[index].map(f64::to_bits),
                    original_stream.f[index].map(f64::to_bits)
                );
            }
        }
    }
    assert!(nonzero_stress > 1e-8);
    let capacity = scratch.capacity();
    for _ in 0..4 {
        grid.step_pressure_x(&mut scratch, boundary);
        assert_eq!(scratch.capacity(), capacity);
    }
}

#[test]
fn pressure_faces_hold_uniform_and_barotropic_hydrostatic_equilibria() {
    // The native Guo/streaming hydrostatic balance is the trapezoidal density
    // relation rho[x+1]/rho[x] = (1+g/(2 cs²))/(1-g/(2 cs²)). This approaches
    // the continuum isothermal law rho(x)=rho(0) exp(g x/cs²) to second order.
    for gravity in [0.0, 4e-5, -4e-5] {
        let mut grid = Grid::uniform(33, 5, 0.9);
        grid.periodic_x = false;
        let ratio: f64 = (1.0 + gravity / (2.0 * CS2)) / (1.0 - gravity / (2.0 * CS2));
        let density = |x: usize| 1.03 * ratio.powi(x as i32);
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                let cell = grid.idx(x, y);
                grid.f[cell] = equilibrium(density(x), 0.0, 0.0);
            }
        }
        grid.rebind_gravity([gravity, 0.0], &mut Vec::new())
            .unwrap();
        let boundary = PressureX2::new(density(0), density(grid.nx - 1));
        let mut scratch = Vec::new();
        for _ in 0..300 {
            grid.step_pressure_x(&mut scratch, boundary);
        }
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                let actual = grid.moments(grid.idx(x, y));
                assert!((actual.rho - density(x)).abs() < 2e-12);
                assert!(
                    actual.u.iter().all(|u| u.abs() < 2e-12),
                    "hydrostatic drift {actual:?}"
                );
                let continuum = 1.03 * (gravity * x as f64 / CS2).exp();
                assert!((actual.rho - continuum).abs() < 1e-11);
            }
        }
    }
}

fn interior_momentum(grid: &Grid) -> f64 {
    (0..grid.ny)
        .flat_map(|y| (1..grid.nx - 1).map(move |x| grid.idx(x, y)))
        .filter(|&cell| grid.flags[cell] == Cell::Fluid)
        .map(|cell| (0..Q).map(|q| EX[q] * grid.f[cell][q]).sum::<f64>())
        .sum()
}

fn open_momentum_flux(grid: &Grid, post: &[[f64; Q]]) -> f64 {
    // Native population transport into the interior control volume. Both
    // diagonal streams are included; periodic y closes transverse transport.
    let mut flux = 0.0;
    for y in 0..grid.ny {
        for (q, &ex) in EX.iter().enumerate() {
            if ex > 0.0 {
                flux += post[grid.idx(0, y)][q] - post[grid.idx(grid.nx - 2, y)][q];
            } else if ex < 0.0 {
                flux += post[grid.idx(1, y)][q] - post[grid.idx(grid.nx - 1, y)][q];
            }
        }
    }
    flux
}

#[test]
fn pressure_head_drives_reversible_flow_and_conserves_wall_momentum() {
    let mut forward = Grid::uniform(21, 10, 0.9);
    forward.periodic_x = false;
    let mut mask = vec![false; forward.f.len()];
    for y in 3..7 {
        let cell = forward.idx(10, y);
        forward.flags[cell] = Cell::Wall;
        mask[cell] = true;
    }
    let mut reverse = forward.clone();
    let head = PressureX2::new(1.00006, 0.99994);
    let reversed_head = PressureX2::new(head.max_density(), head.min_density());
    let (mut scratch, mut reverse_scratch) = (Vec::new(), Vec::new());
    for _ in 0..1600 {
        forward.step_pressure_x(&mut scratch, head);
        reverse.step_pressure_x(&mut reverse_scratch, reversed_head);
    }
    let before = interior_momentum(&forward);
    assert!(
        before > 1e-3,
        "a pressure head must generate flow: {before}"
    );
    forward.collide_into(&mut scratch);
    let pressure_flux = open_momentum_flux(&forward, &scratch);
    let receipt = forward.step_pressure_x_with_wall_momentum(&mut scratch, head, &mask);
    reverse.step_pressure_x(&mut reverse_scratch, reversed_head);
    let after = interior_momentum(&forward);
    assert!(receipt.wall_impulse[0] > 0.0);
    assert!((after - before - pressure_flux + receipt.wall_impulse[0]).abs() < 2e-12);
    assert!((after + interior_momentum(&reverse)).abs() < 2e-11);
    for y in 0..forward.ny {
        for x in 0..forward.nx {
            if forward.flags[forward.idx(x, y)] != Cell::Fluid {
                continue;
            }
            let a = forward.moments(forward.idx(x, y));
            let b = reverse.moments(reverse.idx(reverse.nx - 1 - x, y));
            assert!((a.rho - b.rho).abs() < 2e-12);
            assert!((a.u[0] + b.u[0]).abs() < 2e-12);
            assert!((a.u[1] - b.u[1]).abs() < 2e-12);
        }
    }
}

#[test]
fn invalid_pressure_and_late_source_refusal_do_not_advance_the_grid() {
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(std::panic::catch_unwind(|| PressureX2::new(invalid, 1.0)).is_err());
        assert!(std::panic::catch_unwind(|| PressureX2::new(1.0, invalid)).is_err());
    }
    let mut grid = Grid::uniform(5, 8, 1.0);
    grid.periodic_x = false;
    grid.g = [f64::MAX * 0.5, 0.0];
    let original = grid.f.clone();
    let mut scratch = Vec::new();
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            grid.step_pressure_x(&mut scratch, PressureX2::new(4.0, 1.0));
        }))
        .is_err()
    );
    assert_eq!(grid.f, original);
    assert!(
        scratch.is_empty(),
        "prescribed-density overflow must precede collision"
    );

    grid.g = [0.0; 2];
    for x in 0..grid.nx {
        let cell = grid.idx(x, 5);
        grid.f[cell] = equilibrium(1.0, 0.35, 0.0);
    }
    let original = grid.f.clone();
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            grid.step_pressure_x(&mut scratch, PressureX2::new(1.001, 0.999));
        }))
        .is_err()
    );
    assert_eq!(grid.f, original);
    assert_eq!(scratch.len(), grid.f.len());
}

fn walled_pressure_channel(nx: usize, ny: usize, tau: f64) -> Grid {
    let mut grid = Grid::uniform(nx, ny, tau);
    grid.periodic_x = false;
    grid.periodic_y = false;
    for x in 0..nx {
        for y in [0, ny - 1] {
            let cell = grid.idx(x, y);
            grid.flags[cell] = Cell::Wall;
        }
    }
    grid
}

#[test]
fn constant_profiles_replay_uniform_pressures_and_native_wall_impulses() {
    for sealed_y in [false, true] {
        let mut uniform = walled_pressure_channel(13, 8, 0.9);
        if !sealed_y {
            uniform.periodic_y = true;
            uniform.flags.fill(Cell::Fluid);
        }
        let mut profiled = uniform.clone();
        let min = vec![1.00008; uniform.ny];
        let max = vec![0.99992; uniform.ny];
        let boundary = PressureX2::new(min[0], max[0]);
        let profile = PressureProfileX2::new(&min, &max);
        assert_eq!(profile.min_densities(), min);
        assert_eq!(profile.max_densities(), max);
        let mask: Vec<_> = uniform
            .flags
            .iter()
            .map(|flag| *flag == Cell::Wall)
            .collect();
        let (mut a_scratch, mut b_scratch) = (Vec::new(), Vec::new());
        let mut capacity = None;
        for _ in 0..16 {
            let mut native = profiled.clone();
            let expected = native.step_with_wall_momentum(&mut Vec::new(), &mask);
            let a = uniform.step_pressure_x_with_wall_momentum(&mut a_scratch, boundary, &mask);
            let b =
                profiled.step_pressure_profile_x_with_wall_momentum(&mut b_scratch, profile, &mask);
            assert_eq!(
                a.wall_impulse.map(f64::to_bits),
                b.wall_impulse.map(f64::to_bits)
            );
            assert_eq!(
                b.wall_impulse.map(f64::to_bits),
                expected.wall_impulse.map(f64::to_bits)
            );
            assert_eq!(b.measured_links, expected.measured_links);
            for (a, b) in uniform.f.iter().zip(&profiled.f) {
                assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
            }
            assert_eq!(b_scratch.len(), profiled.f.len());
            if let Some(capacity) = capacity {
                assert_eq!(b_scratch.capacity(), capacity);
            } else {
                capacity = Some(b_scratch.capacity());
            }
        }
    }
}

#[test]
fn transverse_pressure_profiles_preserve_discrete_hydrostatics() {
    // At physical rest, Guo half-source populations satisfy streaming when
    // rho[y+1] / rho[y] = (1 + g/(2 cs²)) / (1 - g/(2 cs²)).
    // Both reservoirs must impose this same transverse pressure field.
    for gravity in [6e-5, -6e-5] {
        let mut grid = walled_pressure_channel(13, 19, 0.9);
        let ratio: f64 = (1.0 + gravity / (2.0 * CS2)) / (1.0 - gravity / (2.0 * CS2));
        let density: Vec<_> = (0..grid.ny).map(|y| 1.01 * ratio.powi(y as i32)).collect();
        for (y, &rho) in density.iter().enumerate() {
            for x in 0..grid.nx {
                let cell = grid.idx(x, y);
                grid.f[cell] = equilibrium(rho, 0.0, 0.0);
            }
        }
        grid.rebind_gravity([0.0, gravity], &mut Vec::new())
            .unwrap();
        let initial = grid.f.clone();
        let boundary = PressureProfileX2::new(&density, &density);
        let mut scratch = Vec::new();
        for _ in 0..400 {
            grid.step_pressure_profile_x(&mut scratch, boundary);
        }
        for (y, &rho) in density.iter().enumerate() {
            for x in 0..grid.nx {
                let cell = grid.idx(x, y);
                if grid.flags[cell] == Cell::Wall {
                    assert_eq!(
                        grid.f[cell].map(f64::to_bits),
                        initial[cell].map(f64::to_bits)
                    );
                } else {
                    let actual = grid.moments(cell);
                    assert!((actual.rho - rho).abs() < 3e-12, "{x},{y}: {actual:?}");
                    assert!(
                        actual.u.iter().all(|u| u.abs() < 3e-12),
                        "{x},{y}: {actual:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn pressure_head_drives_a_no_slip_channel_with_the_analytic_viscosity() {
    // Steady incompressible plane Poiseuille: u(y) = -dp/dx y(H-y)/(2 rho nu).
    // Low head limits compressibility; TRT 3/16 fixes halfway bounce-back walls.
    let mut grid = walled_pressure_channel(33, 10, 0.9);
    grid.collision = CollisionModel2::Trt { magic: 3.0 / 16.0 };
    let min = vec![1.00006; grid.ny];
    let max = vec![0.99994; grid.ny];
    let profile = PressureProfileX2::new(&min, &max);
    let mut scratch = Vec::new();
    for _ in 0..4000 {
        grid.step_pressure_profile_x(&mut scratch, profile);
    }
    let height = (grid.ny - 2) as f64;
    let pressure_gradient = CS2 * (max[0] - min[0]) / (grid.nx - 1) as f64;
    let viscosity = CS2 * (0.9 - 0.5);
    for y in 1..grid.ny - 1 {
        let actual = grid.moments(grid.idx(grid.nx / 2, y));
        let position = y as f64 - 0.5;
        let expected =
            -pressure_gradient * position * (height - position) / (2.0 * actual.rho * viscosity);
        assert!(
            (actual.u[0] / expected - 1.0).abs() < 0.025,
            "row {y}: expected {expected}, got {:?}",
            actual.u
        );
        assert!(actual.u[1].abs() < 2e-8, "transverse drift: {:?}", actual.u);
    }
}

#[test]
fn spatial_pressures_keep_the_local_forcing_frame_and_donor_stress() {
    let mut grid = walled_pressure_channel(9, 8, 0.83);
    let min: Vec<_> = (0..grid.ny).map(|y| 1.002 + 0.0001 * y as f64).collect();
    let max: Vec<_> = (0..grid.ny).map(|y| 0.999 + 0.0002 * y as f64).collect();
    let mut forces = vec![[0.0; 2]; grid.f.len()];
    for y in 1..grid.ny - 1 {
        for x in 0..grid.nx {
            let cell = grid.idx(x, y);
            grid.f[cell] = equilibrium(1.0, 0.01, -0.002);
            grid.f[cell][0] -= 0.0004;
            grid.f[cell][1] += 0.0002;
            grid.f[cell][3] += 0.0002;
            forces[cell] = [1e-5 * x as f64, -7e-6 * y as f64];
        }
    }
    grid.rebind_forcing([2e-5, -3e-5], &forces, &mut Vec::new())
        .unwrap();
    let profile = PressureProfileX2::new(&min, &max);
    let mut scratch = Vec::new();
    for _ in 0..4 {
        grid.step_pressure_profile_x(&mut scratch, profile);
        for y in 1..grid.ny - 1 {
            for (face_x, donor_x, target) in [(0, 1, min[y]), (grid.nx - 1, grid.nx - 2, max[y])] {
                let face = grid.idx(face_x, y);
                let donor = grid.idx(donor_x, y);
                let actual = grid.moments(face);
                let expected = grid.moments(donor);
                assert!((actual.rho - target).abs() < 2e-12);
                for d in 0..2 {
                    assert!((actual.u[d] - expected.u[d]).abs() < 2e-12);
                    for e in 0..2 {
                        assert!(
                            (physical_stress(&grid, face)[d][e]
                                - physical_stress(&grid, donor)[d][e])
                                .abs()
                                < 2e-12
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn pressure_profiles_and_transverse_wall_admission_are_atomic() {
    assert!(std::panic::catch_unwind(|| PressureProfileX2::new(&[], &[])).is_err());
    assert!(std::panic::catch_unwind(|| PressureProfileX2::new(&[1.0], &[1.0, 1.0])).is_err());
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            std::panic::catch_unwind(|| {
                let _ = PressureProfileX2::new(&[invalid], &[1.0]);
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                let _ = PressureProfileX2::new(&[1.0], &[invalid]);
            })
            .is_err()
        );
    }
    let original = walled_pressure_channel(7, 8, 0.9);
    for failure in 0..5 {
        let mut grid = original.clone();
        let mut min = vec![1.001; grid.ny];
        let mut max = vec![0.999; grid.ny];
        match failure {
            0 => {
                min.pop();
                max.pop();
            }
            1 => {
                let cell = grid.idx(3, 0);
                grid.flags[cell] = Cell::Fluid;
            }
            2 => {
                let cell = grid.idx(0, 4);
                grid.flags[cell] = Cell::Wall;
            }
            3 => {
                grid.flags.fill(Cell::Wall);
            }
            _ => {
                for x in 0..grid.nx {
                    let cell = grid.idx(x, grid.ny - 2);
                    grid.f[cell] = equilibrium(1.0, 0.35, 0.0);
                }
            }
        }
        let populations = grid.f.clone();
        let mut scratch = Vec::new();
        let profile = PressureProfileX2::new(&min, &max);
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                grid.step_pressure_profile_x(&mut scratch, profile);
            }))
            .is_err()
        );
        assert_eq!(grid.f, populations);
        if failure < 4 {
            assert!(
                scratch.is_empty(),
                "topology/profile refusal must precede collision"
            );
        } else {
            assert_eq!(scratch.len(), grid.f.len());
        }
    }
}
