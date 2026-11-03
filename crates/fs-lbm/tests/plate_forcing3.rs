//! Physical momentum and stress through held-load changes in dense D3Q19 flow.
use fs_lbm::d3q19::{
    CollisionModel3, E3, ForceRebindError3, PlatesGrid3, Q3, W3, collide_cell3, equilibrium3,
    plate_channel_flow3,
};
use fs_lbm::rheology::Rheology;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn physical_moments(populations: &[f64; Q3], force: [f64; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let rho = populations.iter().sum::<f64>();
    let momentum = core::array::from_fn(|a| {
        populations
            .iter()
            .zip(directions)
            .map(|(f, e)| f * e[a])
            .sum::<f64>()
            + 0.5 * force[a]
    });
    let velocity: [f64; 3] = momentum.map(|p| p / rho);
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            populations
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
        })
    });
    (rho, momentum, stress)
}

fn prepared_grid() -> PlatesGrid3 {
    let mut grid = PlatesGrid3::plates(3, 2, 4, 0.9);
    let stress = [
        [0.0015, -0.0007, 0.0003],
        [-0.0007, -0.0004, 0.0005],
        [0.0003, 0.0005, 0.0002],
    ];
    for cell in 0..grid.f.len() {
        let rho = 0.6 + 0.03 * cell as f64;
        let velocity = [0.02 + 0.0001 * cell as f64, -0.01, 0.015];
        let mut populations = equilibrium3(rho, velocity);
        for (q, (x, y, z)) in E3.into_iter().enumerate() {
            let e = [f64::from(x), f64::from(y), f64::from(z)];
            let contraction = (0..3)
                .map(|a| {
                    (0..3)
                        .map(|b| {
                            stress[a][b] * (e[a] * e[b] - if a == b { 1.0 / 3.0 } else { 0.0 })
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            // Independent Hermite second moment: this cannot be retained by
            // replacing an evolved cell with its equilibrium distribution.
            populations[q] += 4.5 * W3[q] * contraction;
        }
        grid.f[cell] = populations;
        grid.tau[cell] = 0.7 + 0.01 * cell as f64;
    }
    grid
}

#[test]
fn force_changes_preserve_nonuniform_density_momentum_and_full_retained_stress() {
    let mut grid = prepared_grid();
    let initial = grid.f.clone();
    let original: Vec<_> = initial
        .iter()
        .map(|f| physical_moments(f, [0.0; 3]))
        .collect();
    let tau = grid.tau.clone();
    let plane = grid.nx * grid.ny;
    for load in [[0.002, -0.003, 0.004], [-0.003, 0.001, -0.002], [0.0; 3]] {
        let forces: Vec<_> = (0..grid.f.len())
            .map(|cell| load.map(|value| value * (1.0 + 0.01 * cell as f64)))
            .collect();
        grid.rebind_forcing(&forces).unwrap();
        assert_eq!(grid.fext, forces);
        assert_eq!(grid.tau, tau);
        for cell in 0..grid.f.len() {
            if cell < plane || cell >= grid.f.len() - plane {
                assert_eq!(grid.f[cell], initial[cell], "wall populations changed");
                continue;
            }
            let next = physical_moments(&grid.f[cell], grid.fext[cell]);
            near(next.0, original[cell].0, 3e-15);
            for a in 0..3 {
                near(next.1[a], original[cell].1[a], 3e-15);
                for b in 0..3 {
                    near(next.2[a][b], original[cell].2[a][b], 3e-15);
                }
            }
        }
        let accepted = grid.f.clone();
        grid.rebind_forcing(&forces).unwrap();
        assert_eq!(
            grid.f, accepted,
            "identical held forcing changed populations"
        );
    }
}

#[test]
fn switched_force_history_applies_only_the_actual_elapsed_impulse() {
    let rho = 1.7;
    let initial_velocity = [0.015, -0.008, 0.012];
    for tau in [0.6, 0.83, 1.7] {
        let mut grid = PlatesGrid3::plates(2, 3, 24, tau);
        grid.f.fill(equilibrium3(rho, initial_velocity));
        // Four lattice steps cannot reach this cell from either plate. The
        // native collide/stream path must therefore match uniform acceleration.
        let cell = grid.idx(1, 1, 12);
        let mut expected = initial_velocity.map(|u| rho * u);
        for force in [
            [0.0004, -0.0003, 0.0002],
            [-0.0004, 0.0002, 0.0001],
            [0.0; 3],
            [0.0002, 0.0001, -0.0003],
        ] {
            let before = physical_moments(&grid.f[cell], grid.fext[cell]);
            grid.rebind_forcing(&vec![force; grid.f.len()]).unwrap();
            let rebound = physical_moments(&grid.f[cell], force);
            near(rebound.0, before.0, 3e-15);
            for a in 0..3 {
                near(rebound.1[a], before.1[a], 3e-15);
                near(rebound.1[a], expected[a], 1e-14);
            }
            let collided =
                collide_cell3(grid.f[cell], CollisionModel3::Bgk { tau }, force).unwrap();
            let after_collision = physical_moments(&collided, force);
            for a in 0..3 {
                near(after_collision.1[a] - rebound.1[a], force[a], 3e-15);
                expected[a] += force[a];
            }
            grid.step();
            let after = physical_moments(&grid.f[cell], grid.fext[cell]);
            near(after.0, rho, 1e-14);
            for a in 0..3 {
                near(after.1[a], expected[a], 1e-14);
            }
        }
    }
}

#[test]
fn malformed_fields_and_a_late_bad_cell_leave_the_accepted_state_unchanged() {
    let mut grid = prepared_grid();
    let populations = grid.f.clone();
    let forces = grid.fext.clone();
    let tau = grid.tau.clone();
    assert!(matches!(
        grid.rebind_forcing(&vec![[0.001; 3]; grid.f.len() - 1]),
        Err(ForceRebindError3::Invalid(_)),
    ));
    assert_eq!(grid.f, populations);
    assert_eq!(grid.fext, forces);

    let mut proposal = vec![[0.001, -0.002, 0.003]; grid.f.len()];
    let last = proposal.len() - 1;
    proposal[last][2] = f64::NAN;
    assert_eq!(
        grid.rebind_forcing(&proposal),
        Err(ForceRebindError3::Unphysical { cell: last })
    );
    assert_eq!(grid.f, populations);
    assert_eq!(grid.fext, forces);
    proposal[last][2] = 0.003;

    let late = grid.f.len() - grid.nx * grid.ny - 1;
    grid.f[late] = [-1.0; Q3];
    let invalid_populations = grid.f.clone();
    assert_eq!(
        grid.rebind_forcing(&proposal),
        Err(ForceRebindError3::Unphysical { cell: late })
    );
    assert_eq!(
        grid.f, invalid_populations,
        "a staged prefix escaped after late refusal"
    );
    assert_eq!(grid.fext, forces);
    assert_eq!(grid.tau, tau);

    grid.f[late] = populations[late];
    grid.rebind_forcing(&proposal).unwrap();
    let mut clean = prepared_grid();
    clean.rebind_forcing(&proposal).unwrap();
    assert_eq!(grid.f, clean.f, "a refused scratch prefix affected retry");
    assert_eq!(grid.fext, clean.fext);
}

#[test]
fn channel_startup_and_density_perturbation_begin_at_physical_rest() {
    for force in [0.0004, -0.0003] {
        for law in [
            Rheology::Newtonian { nu: 0.07 },
            Rheology::Carreau {
                nu0: 0.2,
                nu_inf: 0.01,
                lambda: 500.0,
                n: 0.6,
            },
        ] {
            let (mut grid, steps, _) = plate_channel_flow3(3, 2, 8, force, law, 0);
            assert_eq!(steps, 0);
            for z in 1..grid.nz - 1 {
                for y in 0..grid.ny {
                    for x in 0..grid.nx {
                        let cell = grid.idx(x, y, z);
                        let (rho, momentum, _) = physical_moments(&grid.f[cell], grid.fext[cell]);
                        near(rho, 1.0, 3e-15);
                        for p in momentum {
                            near(p, 0.0, 3e-15);
                        }
                    }
                }
            }
            let (one_step, steps, _) = plate_channel_flow3(3, 2, 8, force, law, 1);
            assert_eq!(steps, 1);
            let cell = one_step.idx(1, 1, 4);
            let (_, momentum, _) = physical_moments(&one_step.f[cell], one_step.fext[cell]);
            near(momentum[0], force, 3e-15);
            near(momentum[1], 0.0, 3e-15);
            near(momentum[2], 0.0, 3e-15);

            let proposed: Vec<_> = (0..grid.f.len())
                .map(|cell| [force, -0.0002 * (1.0 + cell as f64 / 100.0), 0.0001])
                .collect();
            grid.rebind_forcing(&proposed).unwrap();
            grid.perturb(0x9137_2456, 0.15);
            assert_eq!(grid.fext, proposed);
            for z in 1..grid.nz - 1 {
                for y in 0..grid.ny {
                    for x in 0..grid.nx {
                        let cell = grid.idx(x, y, z);
                        let (rho, momentum, stress) =
                            physical_moments(&grid.f[cell], grid.fext[cell]);
                        assert!((0.85..1.15).contains(&rho));
                        for a in 0..3 {
                            near(momentum[a], 0.0, 3e-15);
                            for b in 0..3 {
                                near(stress[a][b], if a == b { rho / 3.0 } else { 0.0 }, 3e-15);
                            }
                        }
                    }
                }
            }
        }
    }
}
