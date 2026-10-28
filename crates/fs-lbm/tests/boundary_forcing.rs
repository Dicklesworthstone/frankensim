//! Held-source changes must not act as impulses, reset viscous stress, or
//! corrupt the finite-inertia wall balance through a half-force frame jump.

use fs_lbm::d3q19::{CollisionModel3, E3, FiniteMassWall3, Q3, W3};
use fs_lbm::{BoundaryGrid3, BoundarySpec3, CS2, Face3, FaceBoundary3, equilibrium3};

fn cells(grid: &BoundaryGrid3) -> Vec<[usize; 3]> {
    let [nx, ny, nz] = grid.dimensions();
    (0..nz)
        .flat_map(|z| (0..ny).flat_map(move |y| (0..nx).map(move |x| [x, y, z])))
        .filter(|&[x, y, z]| !grid.is_solid(x, y, z))
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

// Independent moment reductions; do not call the production velocity or
// source transform to construct the conservation/stress oracle.
fn moments(grid: &BoundaryGrid3, [x, y, z]: [usize; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let f = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let rho = f.iter().sum::<f64>();
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let momentum = core::array::from_fn(|a| {
        f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]
    });
    let u = momentum.map(|j| j / rho);
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
    (rho, momentum, stress)
}

fn fluid_momentum(grid: &BoundaryGrid3) -> [f64; 3] {
    cells(grid).into_iter().fold([0.0; 3], |sum, cell| {
        let momentum = moments(grid, cell).1;
        core::array::from_fn(|a| sum[a] + momentum[a])
    })
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    cells(grid)
        .into_iter()
        .map(|[x, y, z]| {
            (
                grid.populations(x, y, z).map(f64::to_bits),
                grid.force_density(x, y, z).map(f64::to_bits),
            )
        })
        .collect()
}

fn uniform_grid(model: CollisionModel3) -> BoundaryGrid3 {
    let mut grid =
        BoundaryGrid3::with_collision_model(4, 4, 4, model, [0.0; 3], BoundarySpec3::periodic());
    let initial = equilibrium3(1.4, [0.011, -0.009, 0.006]);
    for [x, y, z] in cells(&grid) {
        grid.set_populations(x, y, z, &initial);
    }
    grid
}

#[test]
fn rebinding_preserves_mass_physical_momentum_and_full_nonequilibrium_stress() {
    let mut grid = BoundaryGrid3::new(8, 4, 4, 0.9, [2e-4, -3e-4, 1e-4], BoundarySpec3::periodic());
    grid.voxelize_sdf(|p| (p[0] - 3.5).abs() - 0.25);
    grid.set_force_field(|[x, y, z]| [x as f64 * 1e-4, y as f64 * -2e-4, z as f64 * 3e-4]);
    for [x, y, z] in cells(&grid) {
        let u = [0.01 + x as f64 * 0.001, -0.02, 0.005];
        let mut f = equilibrium3(1.2 + y as f64 * 0.01, u);
        // Nonzero independent normal and shear stresses, at zero mass and
        // first moment. An equilibrium reset would erase this state.
        f[1] += 0.002;
        f[2] += 0.002;
        f[0] -= 0.004;
        f[7] += 0.001;
        f[8] += 0.001;
        f[9] -= 0.001;
        f[10] -= 0.001;
        let force = grid.force_density(x, y, z);
        for q in 0..Q3 {
            let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
            let eu = e.iter().zip(u).map(|(e, u)| e * u).sum::<f64>();
            let source = (0..3)
                .map(|a| (3.0 * (e[a] - u[a]) + 9.0 * eu * e[a]) * force[a])
                .sum::<f64>();
            f[q] -= 0.5 * W3[q] * source;
        }
        grid.set_populations(x, y, z, &f);
    }
    let before: Vec<_> = cells(&grid)
        .into_iter()
        .map(|cell| moments(&grid, cell))
        .collect();
    for scale in [1.0, -0.7, 0.0] {
        let field = |[x, y, z]: [usize; 3]| {
            assert_ne!(x, 3, "solid cells must not be sampled");
            [
                scale * (y as f64 + 1.0) * 0.001,
                scale * (z as f64 - 1.0) * -0.002,
                scale * (x as f64 + 1.0) * 0.0005,
            ]
        };
        grid.rebind_force_field(field);
        for (cell, expected) in cells(&grid).into_iter().zip(&before) {
            let actual = moments(&grid, cell);
            close(actual.0, expected.0, 3e-15);
            for a in 0..3 {
                close(actual.1[a], expected.1[a], 3e-15);
                for b in 0..3 {
                    close(actual.2[a][b], expected.2[a][b], 3e-15);
                }
            }
        }
        let bits = snapshot(&grid);
        grid.rebind_force_field(field);
        assert_eq!(
            snapshot(&grid),
            bits,
            "unchanged held force must preserve population bits"
        );
    }
}

#[test]
fn switched_reversed_and_disabled_sources_apply_only_elapsed_step_impulses() {
    for model in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 0.9,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 0.9,
            fourth_order_rate: 1.1,
        },
    ] {
        let mut grid = uniform_grid(model);
        let mut expected = fluid_momentum(&grid);
        for force in [[1e-4, -2e-4, 3e-4], [-2e-4, 1e-4, -1e-4], [0.0; 3]] {
            grid.rebind_force_field(|_| force);
            for a in 0..3 {
                close(fluid_momentum(&grid)[a], expected[a], 2e-13);
            }
            for _ in 0..3 {
                grid.step();
                for a in 0..3 {
                    expected[a] += 64.0 * force[a];
                    close(fluid_momentum(&grid)[a], expected[a], 2e-12);
                }
            }
        }
    }
}

#[test]
fn rebinding_and_finite_mass_walls_close_the_external_impulse_and_wall_work_balance() {
    for normal in 0..3 {
        let tangent = [(normal + 1) % 3, (normal + 2) % 3];
        let mut faces = [FaceBoundary3::Periodic; 6];
        for side in 0..2 {
            let mut velocity = [0.0; 3];
            velocity[tangent[0]] = if side == 0 { 0.004 } else { -0.003 };
            faces[2 * normal + side] = FaceBoundary3::Wall { velocity };
        }
        let mut grid = BoundaryGrid3::new(4, 4, 4, 0.8, [0.0; 3], BoundarySpec3::new(faces));
        let walls = [
            FiniteMassWall3 {
                mass: 0.1,
                external_impulse: [2e-5, -3e-5, 4e-5],
            },
            FiniteMassWall3 {
                mass: 1.7,
                external_impulse: [-1e-5, 2e-5, -3e-5],
            },
        ];
        let total = |grid: &BoundaryGrid3| {
            let mut p = fluid_momentum(grid);
            for side in 0..2 {
                let FaceBoundary3::Wall { velocity } =
                    grid.boundaries().face(Face3::ALL[2 * normal + side])
                else {
                    unreachable!()
                };
                for a in 0..3 {
                    p[a] += walls[side].mass * velocity[a];
                }
            }
            p
        };
        for scale in [1.0, -0.5, 0.0] {
            let force = [scale * 1e-5, scale * -2e-5, scale * 3e-5];
            let before = total(&grid);
            grid.rebind_force_field(|_| force);
            for a in 0..3 {
                close(total(&grid)[a], before[a], 2e-13);
            }
            for _ in 0..3 {
                let before = total(&grid);
                let receipts = grid.step_with_finite_mass_walls(walls, [2.0; 3]).unwrap();
                let after = total(&grid);
                for a in 0..3 {
                    let external = 64.0 * force[a]
                        + receipts
                            .iter()
                            .map(|r| r.external_impulse[a] + r.constraint_impulse[a])
                            .sum::<f64>();
                    close(after[a] - before[a], external, 3e-13);
                }
                for receipt in receipts {
                    close(
                        receipt.kinetic_energy_change,
                        receipt.hydrodynamic.wall_work + receipt.external_work,
                        2e-16,
                    );
                }
            }
        }
    }
}

#[test]
fn refused_source_or_late_frame_transform_preserves_force_and_populations_atomically() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let model = CollisionModel3::Bgk { tau: 0.8 };
    let mut grid = uniform_grid(model);
    let mut control = uniform_grid(model);
    grid.rebind_force_field(|_| [1e-4, -2e-4, 3e-4]);
    control.rebind_force_field(|_| [1e-4, -2e-4, 3e-4]);
    grid.step();
    control.step();
    let before = snapshot(&grid);
    for failure in 0..3 {
        assert!(
            catch_unwind(AssertUnwindSafe(|| grid.rebind_force_field(|cell| {
                if cell == [3; 3] {
                    match failure {
                        0 => [f64::NAN, 0.0, 0.0],
                        1 => panic!("late source callback refusal"),
                        _ => [f64::MAX, 0.0, 0.0], // finite force, overflowing source frame
                    }
                } else {
                    [-1e-4, 3e-4, -2e-4]
                }
            })))
            .is_err()
        );
        assert_eq!(snapshot(&grid), before);
    }
    // A refused frame may have touched collision scratch. The next actual
    // collision must still reproduce the untouched control bit for bit.
    grid.step();
    control.step();
    assert_eq!(snapshot(&grid), snapshot(&control));

    let mut invalid = uniform_grid(model);
    invalid.set_populations(3, 3, 3, &[0.0; Q3]);
    let before = snapshot(&invalid);
    assert!(
        catch_unwind(AssertUnwindSafe(
            || invalid.rebind_force_field(|_| [1e-4; 3])
        ))
        .is_err()
    );
    assert_eq!(snapshot(&invalid), before);

    let mut overflow = BoundaryGrid3::new(
        4,
        4,
        4,
        0.8,
        [f64::MAX, 0.0, 0.0],
        BoundarySpec3::periodic(),
    );
    let before = snapshot(&overflow);
    assert!(
        catch_unwind(AssertUnwindSafe(|| overflow.rebind_force_field(|_| [
            f64::MAX,
            0.0,
            0.0
        ])))
        .is_err()
    );
    assert_eq!(snapshot(&overflow), before);
}

#[test]
fn initial_force_frame_semantics_and_open_face_refusal_remain_explicit() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.8, [2e-4, 0.0, 0.0], BoundarySpec3::periodic());
    let populations = grid.populations(1, 1, 1).map(f64::to_bits);
    grid.set_force_field(|_| [0.0, 4e-4, 0.0]);
    assert_eq!(grid.populations(1, 1, 1).map(f64::to_bits), populations);
    close(grid.velocity(1, 1, 1)[0], 1e-4, 1e-16);
    close(grid.velocity(1, 1, 1)[1], 2e-4, 1e-16);

    let mut open = BoundaryGrid3::new(
        4,
        4,
        4,
        0.8,
        [0.0; 3],
        BoundarySpec3::pressure_duct_z(1.0, 1.0),
    );
    let before = snapshot(&open);
    assert!(catch_unwind(AssertUnwindSafe(|| open.rebind_force_field(|_| [0.0; 3]))).is_err());
    assert_eq!(snapshot(&open), before);
}
