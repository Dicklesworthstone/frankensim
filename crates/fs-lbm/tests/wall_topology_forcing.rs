//! Fresh geometry must interpolate physical stress in the held Guo frame.

use fs_lbm::core2::{Cell, Grid};
use fs_lbm::{CS2, Q, equilibrium};

const E: [[f64; 2]; Q] = [
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
const W: [f64; Q] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

#[derive(Clone, Copy, Default)]
struct Physical {
    mass: f64,
    momentum: [f64; 2],
    second: [[f64; 2]; 2],
}

// Independent zeroth/first/second moment identities, without the production
// force-frame transform or Grid::moments as the reference.
fn physical(populations: &[f64; Q], force: [f64; 2]) -> Physical {
    let mass = populations.iter().sum::<f64>();
    let momentum = core::array::from_fn(|a| {
        populations
            .iter()
            .zip(E)
            .map(|(f, e)| f * e[a])
            .sum::<f64>()
            + 0.5 * force[a]
    });
    let velocity = momentum.map(|p| p / mass);
    let second = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            populations
                .iter()
                .zip(E)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
        })
    });
    Physical {
        mass,
        momentum,
        second,
    }
}

fn force(grid: &Grid, cell: usize) -> [f64; 2] {
    let mass = grid.f[cell].iter().sum::<f64>();
    core::array::from_fn(|a| grid.g[a].mul_add(mass, grid.fext[cell][a]))
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-14,
        "actual={actual:e}, expected={expected:e}"
    );
}

fn close_physical(actual: Physical, expected: Physical) {
    close(actual.mass, expected.mass);
    for a in 0..2 {
        close(actual.momentum[a], expected.momentum[a]);
        for b in 0..2 {
            close(actual.second[a][b], expected.second[a][b]);
        }
    }
}

fn fixture(forced: bool) -> (Grid, Vec<bool>, usize) {
    let mut grid = Grid::uniform(5, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = if forced { [0.004, -0.003] } else { [0.0; 2] };
    for cell in 0..grid.f.len() {
        let n = cell as f64;
        let rho = 0.9 + 0.04 * n;
        let u = [0.005 * (n - 3.0), -0.003 * (n - 4.0)];
        let mut h = equilibrium(rho, u[0], u[1]);
        // Anisotropic normal stress and shear are deliberately not equilibrium.
        h[1] += 0.0002 * n;
        h[3] += 0.0002 * n;
        h[0] -= 0.0004 * n;
        h[5] += 0.0001 * n;
        h[7] += 0.0001 * n;
        h[6] -= 0.0001 * n;
        h[8] -= 0.0001 * n;
        grid.tau[cell] = 0.8 + 0.005 * n;
        if forced {
            grid.fext[cell] = [0.002 * (n - 2.0), 0.001 * (7.0 - n)];
            let f: [f64; 2] = core::array::from_fn(|a| rho * grid.g[a] + grid.fext[cell][a]);
            for q in 0..Q {
                let eu = E[q][0] * u[0] + E[q][1] * u[1];
                let ef = E[q][0] * f[0] + E[q][1] * f[1];
                let uf = u[0] * f[0] + u[1] * f[1];
                h[q] -= 0.5 * W[q] * (ef / CS2 + eu * ef / (CS2 * CS2) - uf / CS2);
            }
        }
        grid.f[cell] = h;
    }
    let fresh = grid.idx(1, 1);
    grid.flags[fresh] = Cell::Wall;
    let mut walls = vec![false; grid.f.len()];
    walls[grid.idx(2, 1)] = true;
    (grid, walls, fresh)
}

fn donors(grid: &Grid, walls: &[bool], fresh: usize) -> Vec<usize> {
    let mut result = Vec::new();
    for q in 1..Q {
        if let Some(cell) = grid.source(fresh % grid.nx, fresh / grid.nx, q) {
            if grid.flags[cell] == Cell::Fluid && !walls[cell] && !result.contains(&cell) {
                result.push(cell);
            }
        }
    }
    result
}

#[test]
fn fresh_cells_average_physical_momentum_and_full_stress_before_native_collision() {
    let (mut grid, walls, fresh) = fixture(true);
    let sources = donors(&grid, &walls, fresh);
    assert_eq!(sources.len(), 7);
    let mut expected = Physical::default();
    for &cell in &sources {
        let value = physical(&grid.f[cell], force(&grid, cell));
        expected.mass += value.mass / sources.len() as f64;
        for a in 0..2 {
            expected.momentum[a] += value.momentum[a] / sources.len() as f64;
            for b in 0..2 {
                expected.second[a][b] += value.second[a][b] / sources.len() as f64;
            }
        }
    }
    let _ = grid.transition_wall_topology(&walls);
    let f = force(&grid, fresh);
    close_physical(physical(&grid.f[fresh], f), expected);

    // The reconstructed stress must drive the next real BGK collision. The
    // physical nonequilibrium stress relaxes at 1-1/tau while momentum gets F.
    let u0 = expected.momentum.map(|p| p / expected.mass);
    let u1: [f64; 2] = core::array::from_fn(|a| u0[a] + f[a] / expected.mass);
    let relaxation = 1.0 - 1.0 / grid.tau[fresh];
    let mut post = Vec::new();
    grid.collide_into(&mut post);
    let after = physical(&post[fresh], f);
    close(after.mass, expected.mass);
    for a in 0..2 {
        close(after.momentum[a], expected.momentum[a] + f[a]);
        for b in 0..2 {
            let pressure = if a == b { CS2 } else { 0.0 };
            let eq0 = expected.mass * (pressure + u0[a] * u0[b]);
            let eq1 = expected.mass * (pressure + u1[a] * u1[b]);
            close(
                after.second[a][b],
                eq1 + relaxation * (expected.second[a][b] - eq0),
            );
        }
    }
}

#[test]
fn topology_initialization_is_independent_of_donor_force_representation() {
    let (mut unforced, walls, fresh) = fixture(false);
    let (mut forced, _, _) = fixture(true);
    let _ = unforced.transition_wall_topology(&walls);
    let _ = forced.transition_wall_topology(&walls);
    close_physical(
        physical(&forced.f[fresh], force(&forced, fresh)),
        physical(&unforced.f[fresh], [0.0; 2]),
    );
}

#[test]
fn unforced_topology_retains_original_population_bits_and_unique_periodic_donors() {
    let mut grid = Grid::uniform(2, 2, 0.8);
    grid.flags[0] = Cell::Wall;
    for cell in 1..4 {
        grid.f[cell] = equilibrium(0.7 + 0.2 * cell as f64, 0.01 * cell as f64, -0.02);
    }
    let walls = vec![false; 4];
    let sources = donors(&grid, &walls, 0);
    assert_eq!(sources.len(), 3);
    let mut original = [0.0; Q];
    for &cell in &sources {
        for q in 0..Q {
            original[q] += grid.f[cell][q];
        }
    }
    for f in &mut original {
        *f *= 1.0 / sources.len() as f64;
    }
    let receipt = grid.transition_wall_topology(&walls);
    assert_eq!(receipt.fresh_donor_samples, 3);
    assert_eq!(grid.f[0].map(f64::to_bits), original.map(f64::to_bits));
    let raw = physical(&original, [0.0; 2]);
    close(receipt.fresh_mass, raw.mass);
    for a in 0..2 {
        close(receipt.fresh_momentum[a], raw.momentum[a]);
    }
}

#[test]
fn unrepresentable_topology_force_refuses_without_publishing() {
    for gravity in [[f64::NAN, 0.0], [f64::MAX, 0.0]] {
        let (mut grid, walls, _) = fixture(false);
        grid.g = gravity;
        let original = grid.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = grid.transition_wall_topology(&walls);
        }));
        assert!(result.is_err());
        assert_eq!(grid.flags, original.flags);
        assert_eq!(grid.f, original.f);
        assert_eq!(grid.tau, original.tau);
        assert_eq!(grid.fext, original.fext);
        assert_eq!(grid.g.map(f64::to_bits), original.g.map(f64::to_bits));
    }
}
