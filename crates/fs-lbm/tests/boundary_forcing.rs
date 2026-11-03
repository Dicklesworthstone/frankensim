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
fn scaled_sources_and_finite_mass_walls_close_external_impulse_and_wall_work_balance() {
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
        grid.rebind_force_field(|_| [1e-5, -2e-5, 3e-5]);
        for scale in [1.0, -0.5, 0.0] {
            let force = [scale * 1e-5, scale * -2e-5, scale * 3e-5];
            let before = total(&grid);
            grid.set_force_field_scale(scale);
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
    // Every cell of this 4-deep z duct is an open-face cell or its first
    // interior neighbour: only the force-free field is admitted.
    let before = snapshot(&open);
    assert!(
        catch_unwind(AssertUnwindSafe(
            || open.rebind_force_field(|_| [0.0, 1e-4, 0.0])
        ))
        .is_err()
    );
    assert_eq!(snapshot(&open), before);
    open.rebind_force_field(|_| [0.0; 3]);
    assert_eq!(snapshot(&open), before);
}

fn raw_first_moment(grid: &BoundaryGrid3, keep: impl Fn([usize; 3]) -> bool) -> [f64; 3] {
    cells(grid)
        .into_iter()
        .filter(|&cell| keep(cell))
        .fold([0.0; 3], |sum, [x, y, z]| {
            let f = grid.populations(x, y, z);
            core::array::from_fn(|a| {
                sum[a]
                    + f.iter()
                        .zip(E3)
                        .map(|(f, (ex, ey, ez))| f * f64::from([ex, ey, ez][a]))
                        .sum::<f64>()
            })
        })
}

#[test]
fn open_faces_admit_a_field_that_vanishes_where_they_reconstruct() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let periodic = FaceBoundary3::Periodic;
    let duct = BoundarySpec3::new([
        FaceBoundary3::Velocity {
            velocity: [0.02, 0.0, 0.0],
        },
        FaceBoundary3::Pressure { density: 1.0 },
        periodic,
        periodic,
        periodic,
        periodic,
    ]);
    let pull = 2e-5;
    let mut grid = BoundaryGrid3::new(12, 4, 4, 0.8, [0.0; 3], duct);
    // Cells 0/11 are the open faces and 1/10 their first-interior cells.
    grid.set_force_field(|[x, _, _]| {
        if (2..=9).contains(&x) {
            [0.0, pull, 0.0]
        } else {
            [0.0; 3]
        }
    });
    for reach in [0, 1, 10, 11] {
        let before = snapshot(&grid);
        assert!(
            catch_unwind(AssertUnwindSafe(|| grid.set_force_field(|[x, _, _]| {
                if x == reach {
                    [0.0, pull, 0.0]
                } else {
                    [0.0; 3]
                }
            })))
            .is_err(),
            "a force read by the open face at x={reach} must be refused"
        );
        assert_eq!(snapshot(&grid), before);
    }
    // From rest, cells 1..=10 exchange only rest populations with the faces
    // in the first step, so they gain exactly one held force per forced cell.
    grid.step();
    let gained = raw_first_moment(&grid, |[x, _, _]| (1..=10).contains(&x));
    close(gained[1], 8.0 * 16.0 * pull, 1e-13);
    close(gained[0], 0.0, 1e-13);
}

#[test]
fn scaled_field_preserves_physical_state_and_applies_one_held_force_per_step() {
    let shape = |[x, y, z]: [usize; 3]| {
        [
            1e-5 * (x as f64 + 1.0),
            -2e-5 * (y as f64 + 1.0),
            3e-6 * (z as f64 + 1.0),
        ]
    };
    for model in [
        CollisionModel3::Bgk { tau: 0.7 },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.6,
            higher_order_rate: 1.0,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.6,
            third_order_rate: 1.0,
            fourth_order_rate: 1.0,
        },
    ] {
        let mut scaled = uniform_grid(model);
        let mut direct = uniform_grid(model);
        for grid in [&mut scaled, &mut direct] {
            grid.rebind_force_field(shape);
            for [x, y, z] in cells(grid) {
                let mut f = grid.populations(x, y, z);
                // Carried normal/shear stress with no added mass or momentum.
                f[1] += 0.001;
                f[2] += 0.001;
                f[0] -= 0.002;
                f[7] += 0.0005;
                f[8] += 0.0005;
                f[9] -= 0.0005;
                f[10] -= 0.0005;
                grid.set_populations(x, y, z, &f);
            }
        }
        let total_shape = cells(&scaled).into_iter().fold([0.0; 3], |sum, cell| {
            let local = shape(cell);
            core::array::from_fn(|a| sum[a] + local[a])
        });
        let start = fluid_momentum(&scaled);
        let mut impulse = [0.0; 3];
        for scale in [1.0, 0.6, -0.3, 0.0, 1.7] {
            let before: Vec<_> = cells(&scaled)
                .into_iter()
                .map(|cell| moments(&scaled, cell))
                .collect();
            scaled.set_force_field_scale(scale);
            for (cell, expected) in cells(&scaled).into_iter().zip(before) {
                let actual = moments(&scaled, cell);
                close(actual.0, expected.0, 3e-15);
                for a in 0..3 {
                    close(actual.1[a], expected.1[a], 3e-15);
                    for b in 0..3 {
                        close(actual.2[a][b], expected.2[a][b], 3e-15);
                    }
                }
                assert_eq!(
                    scaled.force_density(cell[0], cell[1], cell[2]),
                    shape(cell).map(|value| scale * value)
                );
            }
            let bits = snapshot(&scaled);
            scaled.set_force_field_scale(scale);
            assert_eq!(
                snapshot(&scaled),
                bits,
                "identical scale must preserve bits"
            );
            direct.rebind_force_field(|cell| shape(cell).map(|value| scale * value));
            assert_eq!(snapshot(&scaled), snapshot(&direct));
            scaled.step();
            direct.step();
            impulse = core::array::from_fn(|a| impulse[a] + scale * total_shape[a]);
            let end = fluid_momentum(&scaled);
            for a in 0..3 {
                close(end[a] - start[a], impulse[a], 3e-13);
            }
        }
        assert_eq!(snapshot(&scaled), snapshot(&direct));
    }
}

#[test]
fn scaled_midpoint_and_endpoint_average_histories_supply_the_prescribed_impulse() {
    let mut grid = uniform_grid(CollisionModel3::Bgk { tau: 0.8 });
    let shape = [2e-4, -1e-4, 3e-5];
    grid.rebind_force_field(|_| shape);
    let mut expected = fluid_momentum(&grid);
    let endpoints = [1.0, 0.6, -0.3, 0.0, 1.7];
    for pair in endpoints.windows(2) {
        // The caller stages the interval average, then restores its accepted
        // endpoint force after stepping. Neither selection supplies an impulse.
        grid.set_force_field_scale(0.5 * (pair[0] + pair[1]));
        for a in 0..3 {
            close(fluid_momentum(&grid)[a], expected[a], 3e-13);
        }
        grid.step();
        grid.set_force_field_scale(pair[1]);
        for a in 0..3 {
            expected[a] += 64.0 * shape[a] * 0.5 * (pair[0] + pair[1]);
            close(fluid_momentum(&grid)[a], expected[a], 3e-13);
        }
    }

    // Independent G1 impulse oracle: integral_0^1 A*t^2 dt = A/3.
    // Refining the caller's source clock reduces midpoint quadrature error
    // by four. A spurious F/2 frame jump would instead be first order.
    let amplitude = 0.003;
    let errors = [4, 8, 16].map(|steps| {
        let mut grid = uniform_grid(CollisionModel3::Bgk { tau: 0.8 });
        grid.rebind_force_field(|_| [amplitude / steps as f64, 0.0, 0.0]);
        grid.set_force_field_scale(0.0);
        let initial = fluid_momentum(&grid)[0];
        for step in 0..steps {
            let midpoint = (step as f64 + 0.5) / steps as f64;
            grid.set_force_field_scale(midpoint * midpoint);
            grid.step();
        }
        let impulse = fluid_momentum(&grid)[0] - initial;
        let error = 64.0 * amplitude / 3.0 - impulse;
        close(
            error,
            64.0 * amplitude / (12.0 * (steps as f64).powi(2)),
            3e-13,
        );
        error
    });
    for pair in errors.windows(2) {
        assert!((pair[0] / pair[1] - 4.0).abs() < 1e-6, "{errors:?}");
    }
}

#[test]
fn invalid_scale_and_late_population_refusal_leave_accepted_state_and_retry_intact() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let model = CollisionModel3::Bgk { tau: 0.8 };
    let mut grid = uniform_grid(model);
    let mut control = uniform_grid(model);
    for grid in [&mut grid, &mut control] {
        grid.rebind_force_field(|_| [1e-4, -2e-4, 3e-4]);
        grid.step();
    }
    let before = snapshot(&grid);
    for scale in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(catch_unwind(AssertUnwindSafe(|| grid.set_force_field_scale(scale))).is_err());
        assert_eq!(snapshot(&grid), before);
    }

    let last = grid.populations(3, 3, 3);
    grid.set_populations(3, 3, 3, &[0.0; Q3]);
    let invalid = snapshot(&grid);
    assert!(catch_unwind(AssertUnwindSafe(|| grid.set_force_field_scale(-0.75))).is_err());
    assert_eq!(
        snapshot(&grid),
        invalid,
        "late refusal must not publish scratch"
    );
    grid.set_populations(3, 3, 3, &last);
    grid.set_force_field_scale(-0.75);
    control.set_force_field_scale(-0.75);
    grid.step();
    control.step();
    assert_eq!(snapshot(&grid), snapshot(&control));

    // The initializer admits finite force data; scaling must separately
    // reject overflow before changing any accepted populations or force.
    let mut overflow = uniform_grid(model);
    overflow.set_force_field(|_| [2.0, 0.0, 0.0]);
    let before = snapshot(&overflow);
    assert!(catch_unwind(AssertUnwindSafe(|| overflow.set_force_field_scale(f64::MAX))).is_err());
    assert_eq!(snapshot(&overflow), before);
}

/// G0/G3/G5: periodic transport cancels, leaving the endpoint trapezoid of
/// the spatially varying force plus the unchanged uniform force. The oracle
/// reduces physical moments independently of production velocity sampling.
#[test]
fn endpoint_force_scale_closes_physical_momentum_and_matches_pooled_bits() {
    use fs_exec::{CancelGate, TilePool};
    let gate = CancelGate::new_clock_free();
    let uniform = [1e-5, -2e-5, 3e-5];
    let shape = |[x, y, z]: [usize; 3]| {
        [
            2e-6 * (x + 1) as f64,
            -3e-6 * (y + 1) as f64,
            1e-6 * (z + 1) as f64,
        ]
    };
    for model in [
        CollisionModel3::Bgk { tau: 0.7 },
        CollisionModel3::Trt {
            tau: 0.8,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.6,
            higher_order_rate: 1.0,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.6,
            third_order_rate: 1.0,
            fourth_order_rate: 1.0,
        },
    ] {
        let make = || {
            let mut grid = BoundaryGrid3::with_collision_model(
                12,
                8,
                8,
                model,
                [0.0; 3],
                BoundarySpec3::periodic(),
            );
            grid.rebind_uniform_force(uniform, &gate).unwrap();
            grid.rebind_force_field(shape);
            grid
        };
        for workers in [1, 3] {
            let mut serial = make();
            let mut pooled = make();
            let mass = serial.total_mass();
            let start = fluid_momentum(&serial);
            let all_cells = cells(&serial);
            let total_shape = all_cells.iter().fold([0.0; 3], |sum, cell| {
                let force = shape(*cell);
                core::array::from_fn(|a| sum[a] + force[a])
            });
            let mut expected = start;
            let mut previous = 1.0;
            let pool = TilePool::for_host(workers, 0xF02CE);
            pool.with_parked_crew_local(|parked| {
                for scale in [0.6, -0.3, 0.0, 1.7] {
                    serial.step_with_endpoint_force_scale(scale, &gate).unwrap();
                    pooled
                        .step_pooled_with_endpoint_force_scale(parked, &gate, scale)
                        .unwrap();
                    expected = core::array::from_fn(|a| {
                        expected[a]
                            + 0.5 * (previous + scale) * total_shape[a]
                            + all_cells.len() as f64 * uniform[a]
                    });
                    let actual = fluid_momentum(&serial);
                    for a in 0..3 {
                        close(actual[a], expected[a], 3e-12);
                    }
                    close(serial.total_mass(), mass, 3e-12);
                    assert_eq!(snapshot(&serial), snapshot(&pooled));
                    previous = scale;
                }
            });
        }
    }
}

/// G1: F(t)=rho*omega*amplitude*cos(omega*t) has the analytic physical
/// velocity increment amplitude*sin(omega*t). Refining endpoint sampling
/// gives second-order convergence at the same quarter-period endpoint.
#[test]
fn endpoint_force_scale_sinusoidal_physical_velocity_converges_quadratically() {
    use fs_exec::CancelGate;
    let gate = CancelGate::new_clock_free();
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut grid = uniform_grid(CollisionModel3::Bgk { tau: 0.9 });
        let amplitude = 0.005;
        let omega = core::f64::consts::FRAC_PI_2 / steps as f64;
        let mass = grid.total_mass();
        let start = fluid_momentum(&grid)[0] / mass;
        grid.rebind_force_field(|_| [1.4 * omega * amplitude, 0.0, 0.0]);
        close(fluid_momentum(&grid)[0] / mass, start, 2e-14);
        for n in 1..=steps {
            grid.step_with_endpoint_force_scale((omega * n as f64).cos(), &gate)
                .unwrap();
        }
        let increment = fluid_momentum(&grid)[0] / mass - start;
        // Independent closed form of the trapezoidal sine sum.
        let discrete = amplitude * (0.5 * omega) / (0.5 * omega).tan();
        close(increment, discrete, 5e-13);
        errors.push((increment - amplitude).abs());
    }
    assert!(errors[0] / errors[1] > 3.99);
    assert!(errors[1] / errors[2] > 3.99);

    // Negative control: replacing a force before collision gives the wrong
    // physical impulse even though its raw population impulse looks correct.
    let mut early = uniform_grid(CollisionModel3::Bgk { tau: 0.9 });
    let force = 1e-4;
    early.rebind_force_field(|_| [force, 0.0, 0.0]);
    let start = fluid_momentum(&early)[0];
    early.set_force_field(|_| [3.0 * force, 0.0, 0.0]);
    early.step();
    let gained = fluid_momentum(&early)[0] - start;
    close(gained, 64.0 * 4.0 * force, 2e-13);
    assert!((gained - 64.0 * 2.0 * force).abs() > 0.01);
}

/// G4: invalid next forcing, a late invalid collision, and cancellation never
/// publish a next force without its population state (or the reverse).
#[test]
fn endpoint_force_scale_refusals_keep_populations_and_force_atomic() {
    use fs_exec::{CancelGate, TilePool};
    use fs_lbm::BoundaryStepError3;
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xF02CE);
    for pooled in [false, true] {
        let advance = |grid: &mut BoundaryGrid3, scale, gate: &CancelGate| {
            if pooled {
                grid.step_pooled_with_endpoint_force_scale(&pool, gate, scale)
            } else {
                grid.step_with_endpoint_force_scale(scale, gate)
            }
        };
        let mut grid = uniform_grid(CollisionModel3::Bgk { tau: 0.9 });
        let before = snapshot(&grid);
        assert!(matches!(
            advance(&mut grid, 1.0, &gate),
            Err(BoundaryStepError3::EndpointForce { .. })
        ));
        assert_eq!(snapshot(&grid), before);
        grid.rebind_force_field(|_| [1e-4, 0.0, 0.0]);
        let before = snapshot(&grid);
        for scale in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                advance(&mut grid, scale, &gate),
                Err(BoundaryStepError3::EndpointForce { .. })
            ));
            assert_eq!(snapshot(&grid), before);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert_eq!(
            advance(&mut grid, 0.3, &cancelled),
            Err(BoundaryStepError3::Cancelled)
        );
        assert_eq!(snapshot(&grid), before);
        grid.set_populations(3, 3, 3, &[0.0; Q3]);
        let before = snapshot(&grid);
        assert!(matches!(
            advance(&mut grid, 0.3, &gate),
            Err(BoundaryStepError3::Collision { .. })
        ));
        assert_eq!(snapshot(&grid), before);

        let mut overflow = uniform_grid(CollisionModel3::Bgk { tau: 0.9 });
        overflow.set_force_field(|_| [2.0, 0.0, 0.0]);
        let before = snapshot(&overflow);
        assert!(matches!(
            advance(&mut overflow, f64::MAX, &gate),
            Err(BoundaryStepError3::EndpointForce { .. })
        ));
        assert_eq!(snapshot(&overflow), before);
    }
}
