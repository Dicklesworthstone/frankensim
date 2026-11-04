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
fn initial_force_frame_semantics_and_open_face_rebinding_remain_explicit() {
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
    // Every cell is on or beside an open face. Rebinding still supplies no
    // impulse and preserves the physical stress even at mixed wall seams.
    let before: Vec<_> = cells(&open)
        .into_iter()
        .map(|p| moments(&open, p))
        .collect();
    for force in [[0.0, 1e-4, 0.0], [0.0; 3]] {
        open.rebind_force_field(|_| force);
        for (p, expected) in cells(&open).into_iter().zip(&before) {
            let actual = moments(&open, p);
            close(actual.0, expected.0, 3e-15);
            for a in 0..3 {
                close(actual.1[a], expected.1[a], 3e-15);
                for b in 0..3 {
                    close(actual.2[a][b], expected.2[a][b], 3e-15);
                }
            }
        }
    }
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
    // From rest, cells 1..=10 exchange only rest populations with the faces
    // in the first step, so they gain exactly one held force per forced cell.
    grid.step();
    let gained = raw_first_moment(&grid, |[x, _, _]| (1..=10).contains(&x));
    close(gained[1], 8.0 * 16.0 * pull, 1e-13);
    close(gained[0], 0.0, 1e-13);
}

fn assert_open_physical_moments(grid: &BoundaryGrid3, axis: usize) {
    let last = grid.dimensions()[axis] - 1;
    for cell in cells(grid) {
        if cell[axis] != 0 && cell[axis] != last {
            continue;
        }
        let upper = cell[axis] == last;
        let mut interior = cell;
        interior[axis] = if upper { last - 1 } else { 1 };
        let (rho, momentum, stress) = moments(grid, cell);
        let (inside_rho, inside_momentum, inside_stress) = moments(grid, interior);
        let face = Face3::ALL[2 * axis + usize::from(upper)];
        let expected = match grid.boundaries().face(face) {
            FaceBoundary3::Velocity { velocity } => {
                close(rho, inside_rho, 2e-14);
                velocity
            }
            FaceBoundary3::Pressure { density } => {
                close(rho, density, 2e-14);
                inside_momentum.map(|p| p / inside_rho)
            }
            _ => unreachable!(),
        };
        for a in 0..3 {
            close(momentum[a] / rho, expected[a], 2e-14);
            for b in 0..3 {
                close(stress[a][b], inside_stress[a][b], 2e-14);
            }
        }
    }
}

/// G0/G3: all six orientations, all velocity components and all nine stress
/// moments use the actual uniform-plus-spatial source at BOTH cells.
#[test]
fn forced_open_faces_prescribe_physical_velocity_density_and_neighbour_stress() {
    use fs_exec::CancelGate;
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for upper in [false, true] {
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * axis + usize::from(upper)] = FaceBoundary3::Velocity {
                velocity: [0.023, -0.017, 0.009],
            };
            faces[2 * axis + usize::from(!upper)] = FaceBoundary3::Pressure { density: 1.07 };
            let mut grid =
                BoundaryGrid3::new(8, 8, 8, 0.8, [2e-4, -3e-4, 1e-4], BoundarySpec3::new(faces));
            grid.set_force_field(|[x, y, z]| {
                [
                    3e-5 * (x + 1) as f64,
                    -2e-5 * (y + 1) as f64,
                    4e-5 * (z + 1) as f64,
                ]
            });
            grid.initialize_uniform_equilibrium(1.1, [0.012, -0.007, 0.006], &gate)
                .unwrap();
            for [x, y, z] in cells(&grid) {
                let mut f = grid.populations(x, y, z);
                // Independent normal and off-diagonal nonequilibrium moments.
                f[1] += 0.001;
                f[2] += 0.001;
                f[0] -= 0.002;
                f[7] += 0.0003;
                f[8] += 0.0003;
                f[9] -= 0.0003;
                f[10] -= 0.0003;
                grid.set_populations(x, y, z, &f);
            }
            grid.step();
            assert_open_physical_moments(&grid, axis);
        }
    }
}

/// G0: forcing changes neither the owner of a wall-crossing link nor a
/// tangentially streamed population at an open/wall intersection.
#[test]
fn forced_open_wall_seams_reconstruct_only_open_owned_links() {
    use fs_exec::CancelGate;
    use fs_lbm::d3q19::{OPP3, collide_cell3};
    let model = CollisionModel3::Bgk { tau: 0.8 };
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [2e-4, -3e-4, 1e-4],
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.03], 1.0),
    );
    grid.set_force_field(|[x, y, z]| [1e-5 * x as f64, -2e-5 * y as f64, 3e-5 * z as f64]);
    grid.initialize_uniform_equilibrium(1.1, [0.01, -0.02, 0.015], &CancelGate::new_clock_free())
        .unwrap();
    let cell = [0, 0, 0];
    let wall = grid.link_mask(0, 0, 0);
    let open = grid.open_link_mask(0, 0, 0);
    assert_ne!(wall, 0);
    assert_ne!(open, 0);
    assert_eq!(wall & open, 0);
    let collision = |[x, y, z]: [usize; 3]| {
        collide_cell3(
            grid.populations(x, y, z),
            model,
            grid.force_density(x, y, z),
        )
        .unwrap()
    };
    let local_post = collision(cell);
    let expected: [f64; Q3] = core::array::from_fn(|q| {
        if wall & (1 << q) != 0 {
            local_post[OPP3[q]]
        } else if open & (1 << q) != 0 {
            0.0 // Checked separately by the physical face-interior oracle.
        } else {
            let e = [E3[q].0, E3[q].1, E3[q].2];
            let source = core::array::from_fn(|a| (cell[a] as isize - e[a] as isize) as usize);
            collision(source)[q]
        }
    });
    grid.step();
    let actual = grid.populations(0, 0, 0);
    for q in 0..Q3 {
        if open & (1 << q) == 0 {
            assert_eq!(actual[q].to_bits(), expected[q].to_bits(), "link {q}");
        }
    }
}

/// G0/G5: held and endpoint forces share physical open traces across serial
/// and pooled stepping. The endpoint check detects an F_old/2 frame leak.
#[test]
fn forced_open_endpoint_history_preserves_targets_and_pooled_bits() {
    use fs_exec::{CancelGate, TilePool};
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0x0F02CE);
    let make = || {
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.02, -0.01, 0.007],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.0 };
        let mut grid =
            BoundaryGrid3::new(8, 8, 8, 0.8, [2e-4, -1e-4, 3e-4], BoundarySpec3::new(faces));
        grid.rebind_force_field(|[x, y, z]| {
            [
                3e-5 * (x + 1) as f64,
                2e-5 * (y + 1) as f64,
                -1e-5 * (z + 1) as f64,
            ]
        });
        grid
    };
    let mut serial = make();
    let mut pooled = make();
    serial.step();
    pooled.step_pooled(&pool, &gate).unwrap();
    assert_eq!(snapshot(&serial), snapshot(&pooled));
    assert_open_physical_moments(&serial, 0);
    for scale in [0.6, -0.4, 0.0, 1.3] {
        serial.step_with_endpoint_force_scale(scale, &gate).unwrap();
        pooled
            .step_pooled_with_endpoint_force_scale(&pool, &gate, scale)
            .unwrap();
        assert_eq!(snapshot(&serial), snapshot(&pooled));
        assert_open_physical_moments(&serial, 0);
    }
}

/// G4: a finite next force can overflow an open-face stress reconstruction.
/// Even a refusal on the second face must preserve both accepted fields.
#[test]
fn forced_open_late_endpoint_overflow_refuses_atomically_and_retries() {
    use fs_exec::{CancelGate, TilePool};
    use fs_lbm::BoundaryStepError3;
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0x0F02CE);
    let make = || {
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.02, 0.0, 0.0],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.0 };
        let mut grid = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::new(faces));
        grid.rebind_force_field(|[x, _, _]| if x >= 6 { [1e-5, 0.0, 0.0] } else { [0.0; 3] });
        grid
    };
    for pooled in [false, true] {
        let advance = |grid: &mut BoundaryGrid3, scale| {
            if pooled {
                grid.step_pooled_with_endpoint_force_scale(&pool, &gate, scale)
            } else {
                grid.step_with_endpoint_force_scale(scale, &gate)
            }
        };
        let mut grid = make();
        let mut control = make();
        let before = snapshot(&grid);
        let error = advance(&mut grid, 1e308).unwrap_err();
        assert!(matches!(error, BoundaryStepError3::Unphysical { cell } if cell[0] >= 6));
        assert_eq!(snapshot(&grid), before);
        advance(&mut grid, 0.5).unwrap();
        advance(&mut control, 0.5).unwrap();
        assert_eq!(snapshot(&grid), snapshot(&control));
        assert_open_physical_moments(&grid, 0);
    }
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

mod endpoint_fields {
    use super::*;
    use fs_exec::{CancelGate, TilePool};
    use fs_lbm::BoundaryStepError3;

    fn uniform(t: f64) -> [f64; 3] {
        [
            2e-5 * (1.0 + t),
            -1e-5 * (1.0 - 2.0 * t),
            3e-5 * (t * t - 0.3),
        ]
    }

    // Independent components evolve with different spatial/time dependences;
    // no scalar scaling of a single initial vector field can represent this.
    fn local([x, y, z]: [usize; 3], t: f64) -> [f64; 3] {
        let (x, y, z) = ((x + 1) as f64, (y + 1) as f64, (z + 1) as f64);
        [
            1e-6 * x * (1.0 - t) + 2e-6 * y * t,
            -1.5e-6 * z * (1.0 + t * t) + 0.5e-6 * x * t,
            1.2e-6 * y * z * t + 0.7e-6 * x,
        ]
    }

    fn prepare(mut grid: BoundaryGrid3) -> BoundaryGrid3 {
        for [x, y, z] in cells(&grid) {
            grid.set_populations(
                x,
                y,
                z,
                &equilibrium3(0.8 + 0.01 * (x + 2 * y + z) as f64, [0.011, -0.007, 0.004]),
            );
        }
        grid.rebind_uniform_force(uniform(0.0), &CancelGate::new_clock_free())
            .unwrap();
        grid.rebind_force_field(|cell| local(cell, 0.0));
        grid
    }

    fn force_sum(grid: &BoundaryGrid3, t: f64) -> [f64; 3] {
        cells(grid).into_iter().fold([0.0; 3], |sum, cell| {
            let local = local(cell, t);
            core::array::from_fn(|a| sum[a] + uniform(t)[a] + local[a])
        })
    }

    #[test]
    fn rotated_spatial_histories_close_impulse_and_match_pooled_bits() {
        let gate = CancelGate::new_clock_free();
        for model in [
            CollisionModel3::Bgk { tau: 0.8 },
            CollisionModel3::Trt {
                tau: 0.9,
                magic: 3.0 / 16.0,
            },
            CollisionModel3::CentralMoment {
                second_order_rate: 1.3,
                higher_order_rate: 0.9,
            },
            CollisionModel3::ReducedCumulant {
                second_order_rate: 1.3,
                third_order_rate: 0.9,
                fourth_order_rate: 1.1,
            },
        ] {
            for workers in [1, 3] {
                let make = || {
                    prepare(BoundaryGrid3::with_collision_model(
                        12,
                        8,
                        8,
                        model,
                        [0.0; 3],
                        BoundarySpec3::periodic(),
                    ))
                };
                let (mut serial, mut pooled) = (make(), make());
                let pool = TilePool::for_host(workers, 0xF02CE);
                let mass = serial.total_mass();
                let mut expected = fluid_momentum(&serial);
                let mut previous = force_sum(&serial, 0.0);
                for t in [-0.4, 0.8, -1.1, 0.2] {
                    let supplied = force_sum(&serial, t);
                    serial
                        .step_with_endpoint_forcing(uniform(t), |cell| local(cell, t), &gate)
                        .unwrap();
                    pooled
                        .step_pooled_with_endpoint_forcing(&pool, &gate, uniform(t), |cell| {
                            local(cell, t)
                        })
                        .unwrap();
                    expected =
                        core::array::from_fn(|a| expected[a] + 0.5 * (previous[a] + supplied[a]));
                    for (actual, expected) in fluid_momentum(&serial).into_iter().zip(expected) {
                        close(actual, expected, 5e-12);
                    }
                    close(serial.total_mass(), mass, 5e-13);
                    assert_eq!(serial.uniform_force_density(), uniform(t));
                    assert_eq!(snapshot(&serial), snapshot(&pooled));
                    previous = supplied;
                }
            }
        }
    }

    #[test]
    fn changing_fields_reconstruct_all_open_faces_in_the_endpoint_frame() {
        let gate = CancelGate::new_clock_free();
        let pool = TilePool::for_host(3, 0xF02CE);
        for axis in 0..3 {
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[2 * axis] = FaceBoundary3::Velocity {
                velocity: [0.02, -0.01, 0.007],
            };
            faces[2 * axis + 1] = FaceBoundary3::Pressure { density: 1.03 };
            let make = || {
                prepare(BoundaryGrid3::new(
                    8,
                    8,
                    8,
                    0.8,
                    [0.0; 3],
                    BoundarySpec3::new(faces),
                ))
            };
            let (mut serial, mut pooled) = (make(), make());
            for t in [-0.4, 0.8, -1.1] {
                serial
                    .step_with_endpoint_forcing(uniform(t), |cell| local(cell, t), &gate)
                    .unwrap();
                pooled
                    .step_pooled_with_endpoint_forcing(&pool, &gate, uniform(t), |cell| {
                        local(cell, t)
                    })
                    .unwrap();
                assert_open_physical_moments(&serial, axis);
                assert_eq!(snapshot(&serial), snapshot(&pooled));
            }
        }
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn angular(grid: &BoundaryGrid3, origin: [f64; 3], force_at: Option<f64>) -> [f64; 3] {
        cells(grid).into_iter().fold([0.0; 3], |sum, cell| {
            let r = core::array::from_fn(|a| cell[a] as f64 + 0.5 - origin[a]);
            let value = force_at.map_or_else(
                || moments(grid, cell).1,
                |t| {
                    let local = local(cell, t);
                    core::array::from_fn(|a| uniform(t)[a] + local[a])
                },
            );
            let torque = cross(r, value);
            core::array::from_fn(|a| sum[a] + torque[a])
        })
    }

    #[test]
    fn changing_force_and_native_wall_loads_close_linear_and_angular_impulse() {
        let gate = CancelGate::new_clock_free();
        let pool = TilePool::for_host(3, 0xF02CE);
        let origin = [-0.3, 0.7, 1.2];
        let make = || {
            let mut grid =
                BoundaryGrid3::new(8, 8, 8, 0.9, [0.0; 3], BoundarySpec3::lid_cavity([0.0; 3]));
            grid.voxelize_sdf(|p| {
                (p[0] - 3.0)
                    .abs()
                    .max((p[1] - 4.0).abs())
                    .max((p[2] - 3.0).abs())
                    - 0.8
            });
            prepare(grid)
        };
        let (mut serial, mut pooled) = (make(), make());
        let links = serial.boundary_links();
        let mass = serial.total_mass();
        let mut previous_t = 0.0;
        for t in [-0.4, 0.8, -1.1] {
            let before = fluid_momentum(&serial);
            let angular_before = angular(&serial, origin, None);
            let force_before = force_sum(&serial, previous_t);
            let force_after = force_sum(&serial, t);
            let torque_before = angular(&serial, origin, Some(previous_t));
            let torque_after = angular(&serial, origin, Some(t));
            let exchange = serial
                .step_with_endpoint_forcing_and_wall_momentum(
                    uniform(t),
                    |cell| local(cell, t),
                    &links,
                    origin,
                    &gate,
                )
                .unwrap();
            let parallel = pooled
                .step_pooled_with_endpoint_forcing_and_wall_momentum(
                    &pool,
                    &gate,
                    uniform(t),
                    |cell| local(cell, t),
                    &links,
                    origin,
                )
                .unwrap();
            let after = fluid_momentum(&serial);
            let angular_after = angular(&serial, origin, None);
            assert!(exchange.measured_links > 0);
            assert_eq!(exchange, parallel);
            assert_eq!(snapshot(&serial), snapshot(&pooled));
            close(serial.total_mass(), mass, 5e-13);
            for a in 0..3 {
                close(
                    after[a] - before[a] + exchange.wall_impulse[a],
                    0.5 * (force_before[a] + force_after[a]),
                    8e-12,
                );
                close(
                    angular_after[a] - angular_before[a] + exchange.wall_angular_impulse[a],
                    0.5 * (torque_before[a] + torque_after[a]),
                    2e-11,
                );
            }
            previous_t = t;
        }
    }

    #[test]
    fn sampled_fields_reduce_to_native_held_and_scalar_history_bits() {
        let gate = CancelGate::new_clock_free();
        // Include open reconstruction: changing only the force metadata after
        // an ordinary step would fail this equality at both faces.
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.02, -0.01, 0.007],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.03 };
        let make = || {
            prepare(BoundaryGrid3::new(
                8,
                4,
                4,
                0.8,
                [0.0; 3],
                BoundarySpec3::new(faces),
            ))
        };
        let (mut sampled, mut reference) = (make(), make());
        for _ in 0..3 {
            sampled
                .step_with_endpoint_forcing(uniform(0.0), |cell| local(cell, 0.0), &gate)
                .unwrap();
            reference.step();
            assert_eq!(snapshot(&sampled), snapshot(&reference));
        }
        for scale in [-0.4, 0.8, 0.0, 1.3] {
            sampled
                .step_with_endpoint_forcing(
                    uniform(0.0),
                    |cell| local(cell, 0.0).map(|f| scale * f),
                    &gate,
                )
                .unwrap();
            reference
                .step_with_endpoint_force_scale(scale, &gate)
                .unwrap();
            assert_eq!(snapshot(&sampled), snapshot(&reference));
        }
        // A newly sampled field is the newly selected spatial shape, including
        // after old scalar scales have cached their previous reference shape.
        sampled
            .step_with_endpoint_forcing(uniform(0.8), |cell| local(cell, 0.8), &gate)
            .unwrap();
        sampled.step_with_endpoint_force_scale(-0.5, &gate).unwrap();
        for [x, y, z] in cells(&sampled) {
            let expected = local([x, y, z], 0.8);
            assert_eq!(
                sampled.force_density(x, y, z),
                core::array::from_fn(|a| uniform(0.8)[a] - 0.5 * expected[a])
            );
        }
    }

    #[test]
    fn sampling_cancellation_and_late_reconstruction_refuse_atomically() {
        let pool = TilePool::for_host(2, 0xF02CE);
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.02, -0.01, 0.007],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.03 };
        let make = || {
            prepare(BoundaryGrid3::new(
                8,
                4,
                4,
                0.8,
                [0.0; 3],
                BoundarySpec3::new(faces),
            ))
        };
        for pooled in [false, true] {
            for case in 0..6 {
                let (mut grid, mut control) = (make(), make());
                let gate = CancelGate::new_clock_free();
                let mut next_uniform = uniform(0.8);
                if case == 0 {
                    next_uniform[1] = f64::NAN;
                }
                if case == 1 {
                    gate.request();
                }
                if case == 5 {
                    grid.set_populations(7, 3, 3, &[0.0; Q3]);
                }
                let before = snapshot(&grid);
                let sample = |cell: [usize; 3]| {
                    if case == 2 && cell == [7, 3, 3] {
                        return [f64::NAN; 3];
                    }
                    if case == 3 && cell == [7, 3, 3] {
                        gate.request();
                    }
                    // Finite sampled data, but the second open face's physical
                    // stress overflows after a successful collision and stream.
                    if case == 4 && cell[0] >= 6 {
                        return [1e303, 0.0, 0.0];
                    }
                    local(cell, 0.8)
                };
                let result = if pooled {
                    grid.step_pooled_with_endpoint_forcing(&pool, &gate, next_uniform, sample)
                } else {
                    grid.step_with_endpoint_forcing(next_uniform, sample, &gate)
                };
                assert!(
                    matches!(
                        result,
                        Err(BoundaryStepError3::Cancelled
                            | BoundaryStepError3::EndpointForce { .. }
                            | BoundaryStepError3::Unphysical { .. }
                            | BoundaryStepError3::Collision { .. })
                    ),
                    "case {case}: {result:?}"
                );
                assert_eq!(snapshot(&grid), before);
                assert_eq!(grid.uniform_force_density(), uniform(0.0));
                if case == 5 {
                    grid.set_populations(7, 3, 3, &control.populations(7, 3, 3));
                }
                let gate = CancelGate::new_clock_free();
                grid.step_with_endpoint_forcing(uniform(0.8), |cell| local(cell, 0.8), &gate)
                    .unwrap();
                control
                    .step_with_endpoint_forcing(uniform(0.8), |cell| local(cell, 0.8), &gate)
                    .unwrap();
                assert_eq!(snapshot(&grid), snapshot(&control));
            }
        }
    }
}

mod linear_resistance {
    use super::*;
    use fs_exec::{CancelGate, TilePool};
    use fs_lbm::BoundaryStepError3;
    use fs_lbm::d3q19::LinearResistance3;

    const AXES: [[f64; 3]; 3] = [
        [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        [-2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
    ];

    fn models() -> [CollisionModel3; 4] {
        [
            CollisionModel3::Bgk { tau: 0.8 },
            CollisionModel3::Trt {
                tau: 0.9,
                magic: 3.0 / 16.0,
            },
            CollisionModel3::CentralMoment {
                second_order_rate: 1.3,
                higher_order_rate: 0.9,
            },
            CollisionModel3::ReducedCumulant {
                second_order_rate: 1.3,
                third_order_rate: 0.9,
                fourth_order_rate: 1.1,
            },
        ]
    }

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a.into_iter().zip(b).map(|(a, b)| a * b).sum()
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    // Independent constitutive oracle, including arbitrary rotation and zero axes.
    fn force(law: LinearResistance3, velocity: [f64; 3]) -> [f64; 3] {
        let relative = core::array::from_fn(|a| law.skeleton_velocity[a] - velocity[a]);
        core::array::from_fn(|a| {
            law.body_force[a]
                + (0..3)
                    .map(|k| law.axes[k][a] * law.principal[k] * dot(law.axes[k], relative))
                    .sum::<f64>()
        })
    }

    fn angular(grid: &BoundaryGrid3, origin: [f64; 3]) -> [f64; 3] {
        cells(grid).into_iter().fold([0.0; 3], |sum, cell| {
            let r = core::array::from_fn(|a| cell[a] as f64 + 0.5 - origin[a]);
            let torque = cross(r, moments(grid, cell).1);
            core::array::from_fn(|a| sum[a] + torque[a])
        })
    }

    #[test]
    fn rotated_uniform_drag_matches_midpoint_momentum_and_work_for_all_collisions() {
        let gate = CancelGate::new_clock_free();
        let rho = 1.4;
        let law = LinearResistance3 {
            principal: [0.04, 0.2, 0.7],
            axes: AXES,
            skeleton_velocity: [0.002, -0.001, 0.0005],
            body_force: [2e-4, -3e-4, 1e-4],
        };
        for model in models() {
            let mut grid = uniform_grid(model);
            let mut expected = [0.011, -0.009, 0.006];
            for _ in 0..5 {
                let previous = expected;
                let relative = core::array::from_fn(|a| previous[a] - law.skeleton_velocity[a]);
                let principal: [f64; 3] = core::array::from_fn(|k| {
                    ((rho - 0.5 * law.principal[k]) * dot(AXES[k], relative)
                        + dot(AXES[k], law.body_force))
                        / (rho + 0.5 * law.principal[k])
                });
                expected = core::array::from_fn(|a| {
                    law.skeleton_velocity[a]
                        + (0..3).map(|k| AXES[k][a] * principal[k]).sum::<f64>()
                });
                let receipt = grid.step_with_linear_resistance(|_| law, &gate).unwrap();
                let actual = grid.velocity(0, 0, 0);
                for a in 0..3 {
                    close(actual[a], expected[a], 2e-14);
                    close(receipt.external_impulse[a], 64.0 * law.body_force[a], 2e-14);
                    close(
                        64.0 * rho * (actual[a] - previous[a]),
                        receipt.external_impulse[a] - receipt.skeleton_impulse[a],
                        2e-13,
                    );
                }
                let kinetic_change = 32.0 * rho * (dot(actual, actual) - dot(previous, previous));
                close(
                    kinetic_change,
                    receipt.external_work + receipt.fluid_drag_work,
                    2e-14,
                );
                close(
                    receipt.fluid_drag_work + receipt.skeleton_work,
                    -receipt.relative_dissipation,
                    2e-14,
                );
                assert!(receipt.relative_dissipation >= 0.0);
                for (actual, expected) in grid
                    .force_density(0, 0, 0)
                    .into_iter()
                    .zip(force(law, actual))
                {
                    close(actual, expected, 2e-14);
                }
            }
        }
    }

    #[test]
    fn strong_resistance_bounds_relative_energy_and_preserves_free_directions() {
        let gate = CancelGate::new_clock_free();
        let rho = 1.4;
        let law = LinearResistance3 {
            principal: [0.0, 4.0 * rho, 100.0 * rho],
            axes: AXES,
            skeleton_velocity: [0.001, -0.002, 0.0005],
            body_force: [0.0; 3],
        };
        let initial_principal = [1e-4, -7e-5, 4e-5];
        let initial = core::array::from_fn(|a| {
            law.skeleton_velocity[a]
                + (0..3)
                    .map(|k| AXES[k][a] * initial_principal[k])
                    .sum::<f64>()
        });
        for model in models() {
            let mut grid = uniform_grid(model);
            grid.initialize_uniform_equilibrium(rho, initial, &gate)
                .unwrap();
            let mut old_energy = 32.0 * rho * dot(initial_principal, initial_principal);
            for step in 1..=12 {
                let receipt = grid.step_with_linear_resistance(|_| law, &gate).unwrap();
                let relative =
                    core::array::from_fn(|a| grid.velocity(0, 0, 0)[a] - law.skeleton_velocity[a]);
                for k in 0..3 {
                    let factor = (rho - 0.5 * law.principal[k]) / (rho + 0.5 * law.principal[k]);
                    close(
                        dot(AXES[k], relative),
                        initial_principal[k] * factor.powi(step),
                        2e-13,
                    );
                }
                let energy = 32.0 * rho * dot(relative, relative);
                assert!(energy <= old_energy + 2e-15);
                close(old_energy - energy, receipt.relative_dissipation, 2e-14);
                if step == 1 {
                    assert!(
                        dot(AXES[1], relative) > 0.0,
                        "stiff midpoint drag may reverse velocity"
                    );
                    assert!(dot(AXES[2], relative) < 0.0);
                }
                old_energy = energy;
            }
        }
    }

    #[test]
    fn resolved_drag_relaxation_converges_at_second_order() {
        let gate = CancelGate::new_clock_free();
        let initial = [0.011, -0.009, 0.006];
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let mut grid = uniform_grid(CollisionModel3::Bgk { tau: 0.8 });
            let law = LinearResistance3 {
                principal: [1.4 * 0.9 / f64::from(steps); 3],
                ..LinearResistance3::default()
            };
            for _ in 0..steps {
                let _ = grid.step_with_linear_resistance(|_| law, &gate).unwrap();
            }
            errors.push((grid.velocity(0, 0, 0)[0] - initial[0] * (-0.9_f64).exp()).abs());
        }
        for pair in errors.windows(2) {
            assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }

    #[test]
    fn transported_endpoint_forces_close_wall_and_skeleton_impulse_and_torque() {
        let gate = CancelGate::new_clock_free();
        let origin = [-0.3, 0.7, 1.2];
        let make = || {
            let mut grid = BoundaryGrid3::new(
                12,
                8,
                8,
                0.9,
                [0.0; 3],
                BoundarySpec3::lid_cavity([0.003, 0.0, 0.0]),
            );
            grid.voxelize_sdf(|p| {
                (p[0] - 3.0)
                    .abs()
                    .max((p[1] - 4.0).abs())
                    .max((p[2] - 3.0).abs())
                    - 0.8
            });
            for [x, y, z] in cells(&grid) {
                grid.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(
                        0.8 + 0.01 * (x + 2 * y + z) as f64,
                        [0.004 + 0.0001 * y as f64, -0.003, 0.001],
                    ),
                );
            }
            grid.rebind_uniform_force([1e-5, -2e-5, 3e-5], &gate)
                .unwrap();
            grid
        };
        let law = |[x, y, z]: [usize; 3], step: usize| LinearResistance3 {
            principal: [
                0.01 + 0.002 * x as f64,
                0.04 + 0.001 * y as f64,
                0.09 + 0.002 * z as f64,
            ],
            axes: AXES,
            skeleton_velocity: [0.0001 * x as f64, -0.001, 0.0002 * z as f64],
            body_force: [
                1e-5 * (y + step + 1) as f64,
                -2e-6 * (z + 1) as f64,
                3e-6 * (x + step) as f64,
            ],
        };
        for workers in [1, 3] {
            let pool = TilePool::for_host(workers, 0xDA2C7);
            let (mut serial, mut pooled) = (make(), make());
            let links = serial.boundary_links();
            for step in 0..4 {
                let before = fluid_momentum(&serial);
                let angular_before = angular(&serial, origin);
                let mass = serial.total_mass();
                let receipt = serial
                    .step_with_linear_resistance_and_wall_momentum(
                        |cell| law(cell, step),
                        &links,
                        origin,
                        &gate,
                    )
                    .unwrap();
                let parallel = pooled
                    .step_pooled_with_linear_resistance_and_wall_momentum(
                        &pool,
                        &gate,
                        |cell| law(cell, step),
                        &links,
                        origin,
                    )
                    .unwrap();
                assert_eq!(snapshot(&serial), snapshot(&pooled));
                assert_eq!(receipt, parallel);
                let after = fluid_momentum(&serial);
                let angular_after = angular(&serial, origin);
                for a in 0..3 {
                    close(
                        after[a] - before[a],
                        receipt.external_impulse[a] - receipt.skeleton_impulse[a]
                            + receipt.walls.fluid_population_impulse[a],
                        2e-11,
                    );
                    close(
                        angular_after[a] - angular_before[a],
                        receipt.external_angular_impulse[a] - receipt.skeleton_angular_impulse[a]
                            + receipt.walls.fluid_population_angular_impulse[a],
                        3e-11,
                    );
                }
                close(
                    serial.total_mass(),
                    mass + receipt.walls.fluid_mass_change,
                    2e-14,
                );
                close(
                    receipt.fluid_drag_work + receipt.skeleton_work,
                    -receipt.relative_dissipation,
                    2e-13,
                );
                for cell @ [x, y, z] in cells(&serial) {
                    let expected = force(law(cell, step), serial.velocity(x, y, z));
                    for (actual, expected) in
                        serial.force_density(x, y, z).into_iter().zip(expected)
                    {
                        close(actual, expected, 2e-14);
                    }
                }
            }
        }
    }

    #[test]
    fn zero_resistance_preserves_native_held_force_population_bits() {
        let gate = CancelGate::new_clock_free();
        let pool = TilePool::for_host(3, 0xDA2C7);
        let body = |[x, y, z]: [usize; 3]| {
            [
                1e-5 * (x + 1) as f64,
                -2e-5 * (y + 1) as f64,
                3e-5 * (z + 1) as f64,
            ]
        };
        for model in models() {
            let make = || {
                let mut grid = uniform_grid(model);
                grid.rebind_force_field(body);
                grid
            };
            let (mut reference, mut serial, mut pooled) = (make(), make(), make());
            for _ in 0..4 {
                reference.step();
                let law = |cell| LinearResistance3 {
                    body_force: body(cell),
                    axes: AXES,
                    ..LinearResistance3::default()
                };
                let receipt = serial.step_with_linear_resistance(law, &gate).unwrap();
                let parallel = pooled
                    .step_pooled_with_linear_resistance(&pool, &gate, law)
                    .unwrap();
                assert_eq!(snapshot(&reference), snapshot(&serial));
                assert_eq!(snapshot(&reference), snapshot(&pooled));
                assert_eq!(receipt, parallel);
                assert_eq!(receipt.skeleton_impulse, [0.0; 3]);
                assert_eq!(receipt.relative_dissipation, 0.0);
            }
        }
    }

    #[test]
    fn tiny_drag_reaction_survives_a_much_larger_independent_body_force() {
        let gate = CancelGate::new_clock_free();
        let mut grid = uniform_grid(models()[0]);
        let initial = grid.velocity(0, 0, 0);
        let momentum_before = fluid_momentum(&grid);
        let law = LinearResistance3 {
            principal: [1e-22, 3e-22, 2e-22],
            axes: AXES,
            skeleton_velocity: [0.001, -0.002, 0.003],
            body_force: [1e-4, -2e-4, 3e-4],
        };
        let receipt = grid.step_with_linear_resistance(|_| law, &gate).unwrap();
        let final_velocity = grid.velocity(0, 0, 0);
        let momentum_after = fluid_momentum(&grid);
        // The constitutive drag is below the total-force ULP in this fixture.
        // Its separate receipt must not conceal a resolved momentum imbalance.
        assert_eq!(grid.force_density(0, 0, 0), law.body_force);
        for a in 0..3 {
            let expected = receipt.external_impulse[a] - receipt.skeleton_impulse[a];
            let scale = momentum_before[a].abs()
                + momentum_after[a].abs()
                + receipt.external_impulse[a].abs()
                + receipt.skeleton_impulse[a].abs();
            assert!(
                (momentum_after[a] - momentum_before[a] - expected).abs()
                    <= 256.0 * f64::EPSILON * scale
            );
        }
        let midpoint = core::array::from_fn(|a| 0.5 * initial[a] + 0.5 * final_velocity[a]);
        let drag_only = LinearResistance3 {
            body_force: [0.0; 3],
            ..law
        };
        for (actual, drag) in receipt
            .skeleton_impulse
            .into_iter()
            .zip(force(drag_only, midpoint))
        {
            let expected = -64.0 * drag;
            assert!(
                (actual - expected).abs() <= 2e-12 * expected.abs(),
                "{actual:e} != {expected:e}"
            );
            assert_ne!(actual, 0.0);
        }
        assert!(receipt.relative_dissipation > 0.0);
        assert!(
            (receipt.fluid_drag_work + receipt.skeleton_work + receipt.relative_dissipation).abs()
                <= 2e-12 * receipt.relative_dissipation
        );
    }

    #[test]
    fn invalid_laws_cancellation_and_late_work_overflow_are_atomic() {
        let pool = TilePool::for_host(3, 0xDA2C7);
        let valid = LinearResistance3 {
            principal: [0.01, 0.02, 0.03],
            axes: AXES,
            ..LinearResistance3::default()
        };
        for pooled in [false, true] {
            for case in 0..7 {
                let (mut grid, mut control) =
                    (uniform_grid(models()[0]), uniform_grid(models()[0]));
                let before = snapshot(&grid);
                let gate = CancelGate::new_clock_free();
                if case == 0 {
                    gate.request();
                }
                let sample = |cell: [usize; 3]| {
                    let mut law = valid;
                    if cell == [3, 3, 3] {
                        match case {
                            1 => law.principal[1] = -0.1,
                            2 => law.axes[2] = law.axes[0],
                            3 => law.skeleton_velocity[2] = f64::NAN,
                            4 => gate.request(),
                            5 => law.body_force = [1e300; 3],
                            _ => {}
                        }
                    }
                    if case == 6 {
                        law.body_force = [0.1, -0.1, 0.1];
                    }
                    law
                };
                let origin = if case == 6 { [1e308; 3] } else { [0.0; 3] };
                let result = if pooled {
                    grid.step_pooled_with_linear_resistance_and_wall_momentum(
                        &pool,
                        &gate,
                        sample,
                        &[],
                        origin,
                    )
                } else {
                    grid.step_with_linear_resistance_and_wall_momentum(sample, &[], origin, &gate)
                };
                assert!(
                    matches!(
                        result,
                        Err(BoundaryStepError3::Cancelled
                            | BoundaryStepError3::EndpointForce { .. }
                            | BoundaryStepError3::Unphysical { .. }
                            | BoundaryStepError3::Collision { .. })
                    ),
                    "case {case}: {result:?}"
                );
                assert_eq!(snapshot(&grid), before);
                let gate = CancelGate::new_clock_free();
                let _ = grid.step_with_linear_resistance(|_| valid, &gate).unwrap();
                let _ = control
                    .step_with_linear_resistance(|_| valid, &gate)
                    .unwrap();
                assert_eq!(snapshot(&grid), snapshot(&control));
            }
        }
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.01, 0.0, 0.0],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.0 };
        let mut grid = BoundaryGrid3::new(8, 4, 4, 0.9, [0.0; 3], BoundarySpec3::new(faces));
        let before = snapshot(&grid);
        assert!(matches!(
            grid.step_with_linear_resistance(
                |_| panic!("must reject open faces before sampling"),
                &CancelGate::new_clock_free()
            ),
            Err(BoundaryStepError3::EndpointForce { .. })
        ));
        assert_eq!(snapshot(&grid), before);
    }

    mod forchheimer {
        use super::*;
        use fs_lbm::d3q19::DarcyForchheimer3;

        fn norm(v: [f64; 3]) -> f64 {
            dot(v, v).sqrt()
        }
        fn nonlinear_force(law: DarcyForchheimer3, velocity: [f64; 3]) -> [f64; 3] {
            let linear = force(law.linear, velocity);
            let v = core::array::from_fn(|a| velocity[a] - law.linear.skeleton_velocity[a]);
            core::array::from_fn(|a| linear[a] - law.quadratic * norm(v) * v[a])
        }

        #[test]
        fn isotropic_quadratic_oracle_closes_work_and_bounds_stiff_relative_energy() {
            let gate = CancelGate::new_clock_free();
            for model in models() {
                for beta in [2.0, 100.0, 8000.0] {
                    let (mut grid, mut rotated) = (uniform_grid(model), uniform_grid(model));
                    let law = DarcyForchheimer3 {
                        linear: LinearResistance3 {
                            principal: [0.08; 3],
                            skeleton_velocity: [0.002, -0.001, 0.0005],
                            ..LinearResistance3::default()
                        },
                        quadratic: beta,
                    };
                    for step in 0..5 {
                        let before = grid.velocity(0, 0, 0);
                        let v0 =
                            core::array::from_fn(|a| before[a] - law.linear.skeleton_velocity[a]);
                        let q = v0.map(|v| (1.4 - 0.04 - 0.5 * beta * norm(v0)) * v);
                        // Direct quadratic formula, independent of production's
                        // scaled conjugate root and material-frame iteration.
                        let speed =
                            ((1.44_f64.powi(2) + 2.0 * beta * norm(q)).sqrt() - 1.44) / beta;
                        let expected = core::array::from_fn::<_, 3, _>(|a| {
                            law.linear.skeleton_velocity[a] + q[a] / (1.44 + 0.5 * beta * speed)
                        });
                        let receipt = grid.step_with_darcy_forchheimer(|_| law, &gate).unwrap();
                        let rotated_law = DarcyForchheimer3 {
                            linear: LinearResistance3 {
                                axes: AXES,
                                ..law.linear
                            },
                            ..law
                        };
                        let same = rotated
                            .step_with_darcy_forchheimer(|_| rotated_law, &gate)
                            .unwrap();
                        assert_eq!(snapshot(&grid), snapshot(&rotated));
                        assert_eq!(receipt, same);
                        let after = grid.velocity(0, 0, 0);
                        let v1 =
                            core::array::from_fn(|a| after[a] - law.linear.skeleton_velocity[a]);
                        for a in 0..3 {
                            close(after[a], expected[a], 3e-13);
                        }
                        let loss = 32.0 * 1.4 * (dot(v0, v0) - dot(v1, v1));
                        assert!(loss >= -2e-14);
                        close(loss, receipt.relative_dissipation, 3e-13);
                        close(
                            receipt.fluid_drag_work + receipt.skeleton_work,
                            -receipt.relative_dissipation,
                            3e-13,
                        );
                        if step == 0 && beta == 8000.0 {
                            assert!(dot(v0, v1) < 0.0);
                        }
                    }
                }
            }
        }

        #[test]
        fn pure_quadratic_decay_is_second_order_and_rotated_anisotropy_is_covariant() {
            let gate = CancelGate::new_clock_free();
            let initial = [0.011, -0.009, 0.006];
            let mut errors = Vec::new();
            for steps in [8, 16, 32] {
                let mut grid = uniform_grid(models()[0]);
                let law = DarcyForchheimer3 {
                    quadratic: 1.4 * 40.0 / f64::from(steps),
                    ..DarcyForchheimer3::default()
                };
                for _ in 0..steps {
                    let _ = grid.step_with_darcy_forchheimer(|_| law, &gate).unwrap();
                }
                errors.push(
                    (grid.velocity(0, 0, 0)[0] - initial[0] / (1.0 + 40.0 * norm(initial))).abs(),
                );
            }
            for pair in errors.windows(2) {
                assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
            }
            let rotate = |v: [f64; 3]| {
                core::array::from_fn::<_, 3, _>(|a| (0..3).map(|k| AXES[k][a] * v[k]).sum())
            };
            let base = DarcyForchheimer3 {
                linear: LinearResistance3 {
                    principal: [0.0, 0.2, 0.7],
                    body_force: [1e-4, -2e-4, 3e-4],
                    skeleton_velocity: [0.001, -0.002, 0.003],
                    ..LinearResistance3::default()
                },
                quadratic: 80.0,
            };
            let rotated = DarcyForchheimer3 {
                linear: LinearResistance3 {
                    axes: AXES,
                    body_force: rotate(base.linear.body_force),
                    skeleton_velocity: rotate(base.linear.skeleton_velocity),
                    ..base.linear
                },
                ..base
            };
            let (mut original, mut transformed) =
                (uniform_grid(models()[0]), uniform_grid(models()[0]));
            transformed
                .initialize_uniform_equilibrium(1.4, rotate(initial), &gate)
                .unwrap();
            for _ in 0..4 {
                let _ = original
                    .step_with_darcy_forchheimer(|_| base, &gate)
                    .unwrap();
                let _ = transformed
                    .step_with_darcy_forchheimer(|_| rotated, &gate)
                    .unwrap();
                for (actual, expected) in transformed
                    .velocity(0, 0, 0)
                    .into_iter()
                    .zip(rotate(original.velocity(0, 0, 0)))
                {
                    close(actual, expected, 3e-14);
                }
            }
        }

        #[test]
        fn anisotropic_transport_closes_wall_skeleton_momentum_and_matches_pool() {
            let gate = CancelGate::new_clock_free();
            let origin = [-0.3, 0.7, 1.2];
            let make = || {
                let mut grid = BoundaryGrid3::new(
                    12,
                    8,
                    8,
                    0.9,
                    [0.0; 3],
                    BoundarySpec3::lid_cavity([0.003, 0.0, 0.0]),
                );
                grid.voxelize_sdf(|p| {
                    (p[0] - 3.0)
                        .abs()
                        .max((p[1] - 4.0).abs())
                        .max((p[2] - 3.0).abs())
                        - 0.8
                });
                for [x, y, z] in cells(&grid) {
                    grid.set_populations(
                        x,
                        y,
                        z,
                        &equilibrium3(0.8 + 0.01 * (x + 2 * y + z) as f64, [0.004, -0.003, 0.001]),
                    );
                }
                grid
            };
            let law = |[x, y, z]: [usize; 3], step: usize| DarcyForchheimer3 {
                linear: LinearResistance3 {
                    principal: [
                        0.02 + 0.001 * x as f64,
                        0.1 + 0.004 * y as f64,
                        0.3 + 0.002 * z as f64,
                    ],
                    axes: AXES,
                    body_force: [
                        1e-5 * (y + step + 1) as f64,
                        -2e-6 * (z + 1) as f64,
                        3e-6 * (x + step) as f64,
                    ],
                    skeleton_velocity: [0.0001 * x as f64, -0.001, 0.0002 * z as f64],
                },
                quadratic: (5 + x + y + z + step) as f64,
            };
            for workers in [1, 3] {
                let pool = TilePool::for_host(workers, 0xF02C4);
                let (mut serial, mut pooled) = (make(), make());
                let links = serial.boundary_links();
                for step in 0..3 {
                    let before = fluid_momentum(&serial);
                    let angular_before = angular(&serial, origin);
                    let receipt = serial
                        .step_with_darcy_forchheimer_and_wall_momentum(
                            |c| law(c, step),
                            &links,
                            origin,
                            &gate,
                        )
                        .unwrap();
                    let parallel = pooled
                        .step_pooled_with_darcy_forchheimer_and_wall_momentum(
                            &pool,
                            &gate,
                            |c| law(c, step),
                            &links,
                            origin,
                        )
                        .unwrap();
                    assert_eq!(snapshot(&serial), snapshot(&pooled));
                    assert_eq!(receipt, parallel);
                    let after = fluid_momentum(&serial);
                    let angular_after = angular(&serial, origin);
                    for a in 0..3 {
                        close(
                            after[a] - before[a],
                            receipt.external_impulse[a] - receipt.skeleton_impulse[a]
                                + receipt.walls.fluid_population_impulse[a],
                            2e-11,
                        );
                        close(
                            angular_after[a] - angular_before[a],
                            receipt.external_angular_impulse[a]
                                - receipt.skeleton_angular_impulse[a]
                                + receipt.walls.fluid_population_angular_impulse[a],
                            3e-11,
                        );
                    }
                    close(
                        receipt.fluid_drag_work + receipt.skeleton_work,
                        -receipt.relative_dissipation,
                        3e-13,
                    );
                    assert!(receipt.relative_dissipation >= 0.0);
                    for c @ [x, y, z] in cells(&serial) {
                        for (actual, expected) in serial
                            .force_density(x, y, z)
                            .into_iter()
                            .zip(nonlinear_force(law(c, step), serial.velocity(x, y, z)))
                        {
                            close(actual, expected, 3e-14);
                        }
                    }
                }
            }
        }

        #[test]
        fn anisotropic_speed_matches_independent_root_at_large_resistance_contrast() {
            let gate = CancelGate::new_clock_free();
            for beta in [0.5, 80.0, 8000.0] {
                let mut grid = uniform_grid(models()[0]);
                grid.initialize_uniform_equilibrium(1.4, [0.01, 1e-8, -0.004], &gate)
                    .unwrap();
                let law = DarcyForchheimer3 {
                    linear: LinearResistance3 {
                        principal: [0.0, 1e6, 0.1],
                        ..LinearResistance3::default()
                    },
                    quadratic: beta,
                };
                for _ in 0..4 {
                    let rho = grid.populations(0, 0, 0).iter().sum::<f64>();
                    let initial = grid.velocity(0, 0, 0);
                    let force = nonlinear_force(law, initial);
                    let q: [f64; 3] = core::array::from_fn(|a| rho * initial[a] + 0.5 * force[a]);
                    let diagonal = law.linear.principal.map(|d| rho + 0.5 * d);
                    // A separate plain bisection starts from the broad
                    // force/rho bound, not production's component-root bracket.
                    let (mut lower, mut upper) = (0.0, norm(q) / rho);
                    for _ in 0..90 {
                        let midpoint = 0.5 * (lower + upper);
                        let velocity =
                            core::array::from_fn(|a| q[a] / (diagonal[a] + 0.5 * beta * midpoint));
                        if norm(velocity) > midpoint {
                            lower = midpoint;
                        } else {
                            upper = midpoint;
                        }
                    }
                    let speed = 0.5 * (lower + upper);
                    let expected: [f64; 3] =
                        core::array::from_fn(|a| q[a] / (diagonal[a] + 0.5 * beta * speed));
                    let _ = grid.step_with_darcy_forchheimer(|_| law, &gate).unwrap();
                    for (actual, expected) in grid.velocity(0, 0, 0).into_iter().zip(expected) {
                        close(actual, expected, 5e-14);
                    }
                }
            }
        }

        #[test]
        fn zero_quadratic_reduces_exactly_and_bad_nonlinear_endpoints_are_atomic() {
            let gate = CancelGate::new_clock_free();
            let pool = TilePool::for_host(3, 0xF02C4);
            let linear = LinearResistance3 {
                principal: [0.01, 0.04, 0.1],
                axes: AXES,
                body_force: [1e-5, -2e-5, 3e-5],
                ..LinearResistance3::default()
            };
            for model in models() {
                let (mut original, mut extension) = (uniform_grid(model), uniform_grid(model));
                for _ in 0..3 {
                    let receipt = original
                        .step_with_linear_resistance(|_| linear, &gate)
                        .unwrap();
                    let reduced = extension
                        .step_pooled_with_darcy_forchheimer(&pool, &gate, |_| linear.into())
                        .unwrap();
                    assert_eq!(snapshot(&original), snapshot(&extension));
                    assert_eq!(receipt, reduced);
                }
            }
            for pooled in [false, true] {
                for case in 0..4 {
                    let (mut grid, mut control) =
                        (uniform_grid(models()[0]), uniform_grid(models()[0]));
                    if case == 3 {
                        grid.initialize_uniform_equilibrium(1e-300, [0.01, -0.02, 0.03], &gate)
                            .unwrap();
                        control
                            .initialize_uniform_equilibrium(1e-300, [0.01, -0.02, 0.03], &gate)
                            .unwrap();
                    }
                    let skeleton = grid.velocity(0, 0, 0);
                    let before = snapshot(&grid);
                    let cancel = CancelGate::new_clock_free();
                    let sample = |cell: [usize; 3]| {
                        let mut law = DarcyForchheimer3 {
                            linear,
                            quadratic: 10.0,
                        };
                        if cell == [3, 3, 3] {
                            match case {
                                0 => law.quadratic = -1.0,
                                1 => law.quadratic = f64::NAN,
                                2 => cancel.request(),
                                _ => {}
                            }
                        }
                        if case == 3 {
                            // Initial relative velocity is zero: collision is
                            // representable, but the endpoint principal speed
                            // q/(rho+D/2) underflows and must not publish.
                            law.linear = LinearResistance3 {
                                principal: [f64::MAX, 0.0, 0.0],
                                body_force: [1e-301, 0.0, 0.0],
                                skeleton_velocity: skeleton,
                                ..LinearResistance3::default()
                            };
                        }
                        law
                    };
                    let result = if pooled {
                        grid.step_pooled_with_darcy_forchheimer(&pool, &cancel, sample)
                    } else {
                        grid.step_with_darcy_forchheimer(sample, &cancel)
                    };
                    assert!(
                        matches!(
                            result,
                            Err(BoundaryStepError3::Cancelled
                                | BoundaryStepError3::EndpointForce { .. }
                                | BoundaryStepError3::Unphysical { .. })
                        ),
                        "case {case}: {result:?}"
                    );
                    assert_eq!(snapshot(&grid), before);
                    let _ = grid
                        .step_with_darcy_forchheimer(|_| DarcyForchheimer3::default(), &gate)
                        .unwrap();
                    let _ = control
                        .step_with_darcy_forchheimer(|_| DarcyForchheimer3::default(), &gate)
                        .unwrap();
                    assert_eq!(snapshot(&grid), snapshot(&control));
                }
            }
        }
    }
}
