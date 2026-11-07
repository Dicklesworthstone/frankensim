//! Analytic resolved-vorticity checks on the original planar D2Q9 owner.

use fs_exec::CancelGate;
use fs_lbm::core2::{Cell, Grid};
use fs_lbm::equilibrium;
use fs_lbm::kinematics::{KinematicsError, KinematicsWorkspace3};

fn workspace(grid: &Grid) -> KinematicsWorkspace3 {
    let dims = [grid.nx, grid.ny, 1];
    KinematicsWorkspace3::new(dims, KinematicsWorkspace3::required_bytes(dims).unwrap()).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tol {tolerance})"
    );
}

fn prescribe(grid: &mut Grid, velocity: impl Fn(usize, usize) -> [f64; 2]) {
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            let i = grid.idx(x, y);
            let u = velocity(x, y);
            let rho = 1.0 + 0.001 * x as f64;
            let acceleration = grid.fext[i].map(|f| f / rho);
            grid.f[i] = equilibrium(
                rho,
                u[0] - 0.5 * (grid.g[0] + acceleration[0]),
                u[1] - 0.5 * (grid.g[1] + acceleration[1]),
            );
        }
    }
}

#[test]
fn rotation_strain_shear_and_dilation_retain_the_full_tensor() {
    let gate = CancelGate::new_clock_free();
    for (g, curl, divergence, q) in [
        ([[0.0, -0.002], [0.002, 0.0]], 0.004, 0.0, 4e-6),
        ([[0.002, 0.001], [0.001, -0.002]], 0.0, 0.0, -5e-6),
        ([[0.0, 0.003], [0.0, 0.0]], -0.003, 0.0, 0.0),
        ([[0.002, 0.0], [0.0, 0.002]], 0.0, 0.004, -4e-6),
    ] {
        for boost in [[0.0; 2], [0.017, -0.009]] {
            let mut grid = Grid::uniform(7, 7, 0.8);
            grid.periodic_x = false;
            grid.periodic_y = false;
            grid.g = [0.001, -0.002];
            for y in 0..grid.ny {
                for x in 0..grid.nx {
                    let i = grid.idx(x, y);
                    grid.fext[i] = [0.0002 * y as f64, 0.0003 * x as f64];
                }
            }
            prescribe(&mut grid, |x, y| {
                [
                    boost[0] + g[0][0] * x as f64 + g[0][1] * y as f64,
                    boost[1] + g[1][0] * x as f64 + g[1][1] * y as f64,
                ]
            });
            let before = grid.f.clone();
            let mut w = workspace(&grid);
            let flow = grid.kinematics(&gate, &mut w).unwrap()[grid.idx(3, 3)].unwrap();
            for (row, expected) in flow.gradient().iter().take(2).zip(g) {
                for (&actual, expected) in row.iter().take(2).zip(expected) {
                    near(actual, expected, 8e-17);
                }
            }
            near(flow.vorticity()[2], curl, 1e-16);
            near(flow.divergence(), divergence, 1e-16);
            near(flow.q_criterion(), q, 1e-18);
            near(flow.enstrophy_density(), 0.5 * curl * curl, 1e-18);
            assert_eq!(flow.vortex_stretching(), [0.0; 3]);
            near(
                flow.to_physical(0.002).unwrap().vorticity()[2],
                curl / 0.002,
                1e-13,
            );
            assert_eq!(grid.f, before);
        }
    }
}

#[test]
fn moving_halfway_walls_recover_both_orientations_and_one_cell_gaps() {
    let gate = CancelGate::new_clock_free();
    for normal in 0..2 {
        let tangent = 1 - normal;
        for height in [1, 7] {
            for (linear, quadratic) in [(0.0, 0.0), (0.002, 0.0), (0.001, 0.0003)] {
                let mut dims = [5; 2];
                dims[normal] = height + 2;
                let mut grid = Grid::uniform(dims[0], dims[1], 0.8);
                grid.periodic_x = normal != 0;
                grid.periodic_y = normal != 1;
                let profile = |s: f64| 0.02 + linear * s + quadratic * s * s;
                let mut walls = vec![[0.0; 2]; grid.f.len()];
                for y in 0..grid.ny {
                    for x in 0..grid.nx {
                        let i = grid.idx(x, y);
                        let p = [x, y][normal];
                        if p == 0 || p == height + 1 {
                            grid.flags[i] = Cell::Wall;
                            walls[i][tangent] = profile(if p == 0 { 0.0 } else { height as f64 });
                        } else {
                            grid.fext[i][tangent] = 0.0005 * (p as f64 - 0.5);
                        }
                    }
                }
                prescribe(&mut grid, |x, y| {
                    let mut u = [0.0; 2];
                    u[tangent] = profile([x, y][normal] as f64 - 0.5);
                    u
                });
                for i in 0..grid.f.len() {
                    if grid.flags[i] == Cell::Wall {
                        grid.f[i].fill(f64::NAN); // Solid populations are not a velocity sample.
                    }
                }
                let mut w = workspace(&grid);
                let field = grid
                    .kinematics_with_moving_walls(&gate, &walls, &mut w)
                    .unwrap();
                for y in 0..grid.ny {
                    for x in 0..grid.nx {
                        let i = grid.idx(x, y);
                        if grid.flags[i] == Cell::Wall {
                            assert!(field[i].is_none());
                            continue;
                        }
                        let s = [x, y][normal] as f64 - 0.5;
                        let derivative = linear + 2.0 * quadratic * s;
                        let flow = field[i].unwrap();
                        near(flow.gradient()[tangent][normal], derivative, 2e-16);
                        near(
                            flow.vorticity()[2],
                            if normal == 0 { derivative } else { -derivative },
                            2e-16,
                        );
                        near(flow.q_criterion(), 0.0, 1e-18);
                        assert_eq!(flow.vortex_stretching(), [0.0; 3]);
                    }
                }
            }
        }
    }
}

#[test]
fn outer_stationary_halfway_faces_recover_quadratic_channel_curl() {
    let mut grid = Grid::uniform(1, 9, 0.8); // One periodic cell is a homogeneous x extrusion.
    grid.periodic_y = false;
    prescribe(&mut grid, |_, y| {
        let s = y as f64 + 0.5;
        [0.001 * s * (9.0 - s), 0.0]
    });
    let mut w = workspace(&grid);
    let field = grid
        .kinematics(&CancelGate::new_clock_free(), &mut w)
        .unwrap();
    for (y, flow) in field.iter().enumerate() {
        near(
            flow.unwrap().vorticity()[2],
            -0.001 * (9.0 - 2.0 * (y as f64 + 0.5)),
            2e-16,
        );
    }
}

fn taylor_green(n: usize) -> Grid {
    let mut grid = Grid::uniform(n, n, 0.8);
    // Constant density avoids injecting an unrelated density perturbation.
    for y in 0..n {
        for x in 0..n {
            let k = core::f64::consts::TAU / n as f64;
            let (sx, cx) = (k * (x as f64 + 0.5)).sin_cos();
            let (sy, cy) = (k * (y as f64 + 0.5)).sin_cos();
            let i = grid.idx(x, y);
            grid.f[i] = equilibrium(1.0, 0.01 * sx * cy, -0.01 * cx * sy);
        }
    }
    grid
}

#[test]
fn periodic_taylor_green_curl_converges_at_second_order() {
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let grid = taylor_green(n);
        let mut w = workspace(&grid);
        let field = grid
            .kinematics(&CancelGate::new_clock_free(), &mut w)
            .unwrap();
        let mut square_error = 0.0;
        let mut circulation = 0.0;
        for y in 0..n {
            for x in 0..n {
                let k = core::f64::consts::TAU;
                let expected = 0.02
                    * k
                    * (k * (x as f64 + 0.5) / n as f64).sin()
                    * (k * (y as f64 + 0.5) / n as f64).sin();
                let flow = field[grid.idx(x, y)]
                    .unwrap()
                    .to_physical(1.0 / n as f64)
                    .unwrap();
                square_error += (flow.vorticity()[2] - expected).powi(2);
                circulation += flow.vorticity()[2] / (n * n) as f64;
                near(flow.divergence(), 0.0, 5e-15);
                assert_eq!(flow.vortex_stretching(), [0.0; 3]);
            }
        }
        near(circulation, 0.0, 1e-16);
        errors.push((square_error / (n * n) as f64).sqrt());
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.9, "curl errors: {errors:?}");
    }
}

#[test]
fn native_viscosity_decays_taylor_green_enstrophy() {
    let mut grid = taylor_green(32);
    let gate = CancelGate::new_clock_free();
    let mut w = workspace(&grid);
    let enstrophy = |field: &[Option<fs_lbm::kinematics::FlowKinematics3>]| {
        field
            .iter()
            .map(|f| f.unwrap().enstrophy_density())
            .sum::<f64>()
    };
    let initial = enstrophy(grid.kinematics(&gate, &mut w).unwrap());
    let mut scratch = Vec::new();
    for _ in 0..64 {
        grid.step(&mut scratch);
    }
    let final_value = enstrophy(grid.kinematics(&gate, &mut w).unwrap());
    let nu = (0.8 - 0.5) / 3.0;
    let k = core::f64::consts::TAU / 32.0;
    let expected_ratio = (-4.0 * nu * k * k * 64.0).exp();
    near(final_value / initial / expected_ratio, 1.0, 0.02);
}

#[test]
fn incomplete_or_unresolved_fields_refuse_without_changing_the_flow() {
    let mut grid = Grid::uniform(5, 5, 0.8);
    let mut w = workspace(&grid);
    let gate = CancelGate::new_clock_free();
    let expected = grid.kinematics(&gate, &mut w).unwrap().to_vec();
    let before = grid.f.clone();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        grid.kinematics(&cancelled, &mut w),
        Err(KinematicsError::Cancelled)
    );
    for flag in [Cell::Gas, Cell::Interface] {
        grid.flags[12] = flag;
        assert!(grid.kinematics(&gate, &mut w).is_err());
    }
    grid.flags[12] = Cell::Fluid;
    let mut invalid_wall = vec![[0.0; 2]; grid.f.len()];
    invalid_wall[0] = [0.01, 0.0];
    assert!(
        grid.kinematics_with_moving_walls(&gate, &invalid_wall, &mut w)
            .is_err()
    );
    assert!(
        grid.kinematics_with_moving_walls(&gate, &[], &mut w)
            .is_err()
    );
    grid.f[12] = [0.0; 9];
    assert_eq!(
        grid.kinematics(&gate, &mut w),
        Err(KinematicsError::InvalidState { cell: [2, 2, 0] })
    );
    grid.f[12] = before[12];
    assert_eq!(grid.kinematics(&gate, &mut w).unwrap(), expected);
    assert_eq!(grid.f, before);
    let tiny = Grid::uniform(2, 5, 0.8);
    assert!(tiny.kinematics(&gate, &mut workspace(&tiny)).is_err());
    assert!(grid.kinematics(&gate, &mut workspace(&tiny)).is_err());
}
