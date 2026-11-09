//! G0/G1/G3/G4/G5: density-coupled acceleration at native reservoir boundaries.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryStepError3, CollisionModel3, E3, Face3, InletFlux3, OPP3, Q3, RigidFrame3,
    collide_cell3,
};
use fs_lbm::{BoundaryGrid3, BoundarySpec3, FaceBoundary3, equilibrium3};

const N: usize = 8;

fn cells() -> impl Iterator<Item = [usize; 3]> {
    (0..N).flat_map(|z| (0..N).flat_map(move |y| (0..N).map(move |x| [x, y, z])))
}

fn index([x, y, z]: [usize; 3]) -> usize {
    (z * N + y) * N + x
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn model(face: Face3) -> CollisionModel3 {
    [
        CollisionModel3::Bgk { tau: 0.9 },
        CollisionModel3::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.2,
            higher_order_rate: 0.9,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.2,
            third_order_rate: 0.9,
            fourth_order_rate: 1.1,
        },
    ][face as usize % 4]
}

fn moments(grid: &BoundaryGrid3, [x, y, z]: [usize; 3]) -> (f64, [f64; 3]) {
    let populations = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let rho = populations.iter().sum::<f64>();
    let momentum = core::array::from_fn(|a| {
        populations
            .iter()
            .zip(E3)
            .map(|(f, e)| f * f64::from([e.0, e.1, e.2][a]))
            .sum::<f64>()
            + 0.5 * force[a]
    });
    (rho, momentum)
}

fn totals(grid: &BoundaryGrid3) -> (f64, [f64; 3], [f64; 3]) {
    let mut result = (0.0, [0.0; 3], [0.0; 3]);
    for cell @ [x, y, z] in cells() {
        let (rho, momentum) = moments(grid, cell);
        let force = grid.force_density(x, y, z);
        result.0 += rho;
        for a in 0..3 {
            result.1[a] += momentum[a];
            result.2[a] += force[a];
        }
    }
    result
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    cells()
        .map(|[x, y, z]| {
            (
                grid.populations(x, y, z).map(f64::to_bits),
                grid.force_density(x, y, z).map(f64::to_bits),
            )
        })
        .collect()
}

#[test]
fn all_face_normals_keep_exact_uniform_endpoint_response_at_unequal_densities() {
    let gate = CancelGate::new_clock_free();
    for face in Face3::ALL {
        for density in [0.4, 2.3] {
            let mut velocity = [0.011, -0.008, 0.006];
            let mut previous = [0.0007, -0.0012, 0.0009];
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[face as usize] = FaceBoundary3::Velocity { velocity };
            faces[(face as usize) ^ 1] = FaceBoundary3::Pressure { density };
            let mut grid = BoundaryGrid3::with_collision_model(
                N,
                N,
                N,
                model(face),
                [0.0; 3],
                BoundarySpec3::new(faces),
            );
            grid.initialize_uniform_equilibrium(density, velocity, &gate)
                .unwrap();
            let before_bind = totals(&grid);
            grid.rebind_body_acceleration(previous, &gate).unwrap();
            let after_bind = totals(&grid);
            close(after_bind.0, before_bind.0, 2e-11);
            for a in 0..3 {
                close(after_bind.1[a], before_bind.1[a], 2e-12);
            }
            for acceleration in [previous, [-0.0003, 0.0004, 0.0008], [0.0; 3]] {
                for a in 0..3 {
                    velocity[a] += 0.5 * (previous[a] + acceleration[a]);
                }
                faces[face as usize] = FaceBoundary3::Velocity { velocity };
                grid.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
                    .unwrap();
                let before = totals(&grid);
                grid.step_with_endpoint_body_acceleration(acceleration, &gate)
                    .unwrap();
                for cell @ [x, y, z] in cells() {
                    let (rho, momentum) = moments(&grid, cell);
                    close(rho, density, 5e-14);
                    for a in 0..3 {
                        close(momentum[a] / rho, velocity[a], 2e-14);
                        assert_eq!(grid.force_density(x, y, z)[a], rho * acceleration[a]);
                    }
                }
                // This uniform reservoir trace matches the accelerated interior:
                // it needs neither extra mass nor extra physical impulse.
                let after = totals(&grid);
                close(after.0, before.0, 3e-10);
                for a in 0..3 {
                    close(
                        after.1[a] - before.1[a],
                        0.5 * (before.2[a] + after.2[a]),
                        3e-11,
                    );
                }
                previous = acceleration;
            }
        }
    }
}

fn velocity_trace(cell: [usize; 3], axis: usize) -> [f64; 3] {
    [0.012 + 0.0001 * cell[(axis + 1) % 3] as f64, -0.006, 0.004]
}

fn pressure_trace(cell: [usize; 3], axis: usize) -> f64 {
    1.02 + 0.001 * cell[(axis + 2) % 3] as f64
}

fn seamed_grid(face: Face3, gate: &CancelGate) -> BoundaryGrid3 {
    let axis = face as usize / 2;
    let opposite = Face3::ALL[(face as usize) ^ 1];
    let mut faces = [FaceBoundary3::stationary_wall(); 6];
    faces[face as usize] = FaceBoundary3::Velocity {
        velocity: velocity_trace([0; 3], axis),
    };
    faces[opposite as usize] = FaceBoundary3::Pressure {
        density: pressure_trace([0; 3], axis),
    };
    let mut grid = BoundaryGrid3::with_collision_model(
        N,
        N,
        N,
        model(face),
        [0.0; 3],
        BoundarySpec3::new(faces),
    );
    for [x, y, z] in cells() {
        let mut populations = equilibrium3(
            0.9 + 0.004 * x as f64 + 0.003 * y as f64 + 0.002 * z as f64,
            [
                0.008 + 0.0003 * y as f64,
                -0.004 + 0.0002 * z as f64,
                0.003 - 0.0001 * x as f64,
            ],
        );
        // A retained physical shear mode makes a stress reset detectable.
        populations[7] += 1e-5;
        populations[8] += 1e-5;
        populations[9] -= 1e-5;
        populations[10] -= 1e-5;
        grid.set_populations(x, y, z, &populations);
    }
    grid.set_velocity_profile(face, InletFlux3::Unconstrained, gate, |cell| {
        Ok(velocity_trace(cell, axis))
    })
    .unwrap();
    grid.set_pressure_profile(opposite, gate, |cell| Ok(pressure_trace(cell, axis)))
        .unwrap();
    grid.rebind_body_acceleration([0.001, -0.0007, 0.0004], gate)
        .unwrap();
    grid
}

// Independently converge the original prescribed-force owner from the SAME
// accepted state. Each trial supplies fixed force density, not the new scalar
// density closure or an extra elapsed reframe/collision on the previous trial.
fn prescribed_force_reference(
    initial: &BoundaryGrid3,
    acceleration: impl Fn([usize; 3]) -> [f64; 3],
    gate: &CancelGate,
) -> BoundaryGrid3 {
    let mut density: Vec<_> = cells().map(|cell| moments(initial, cell).0).collect();
    for _ in 0..20 {
        let mut next = initial.clone();
        next.step_with_endpoint_forcing(
            [0.0; 3],
            |cell| acceleration(cell).map(|a| a * density[index(cell)]),
            gate,
        )
        .unwrap();
        let updated: Vec<_> = cells().map(|cell| moments(&next, cell).0).collect();
        let error = density
            .iter()
            .zip(&updated)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        if error < 8e-15 {
            return next;
        }
        density = updated;
    }
    panic!("prescribed-force reference failed to converge");
}

fn streamed_before_reconstruction(
    grid: &BoundaryGrid3,
    post: &[[f64; Q3]],
    cell @ [x, y, z]: [usize; 3],
) -> [f64; Q3] {
    let mask = grid.link_mask(x, y, z) | grid.open_link_mask(x, y, z);
    core::array::from_fn(|q| {
        if mask & (1 << q) != 0 {
            post[index(cell)][OPP3[q]]
        } else {
            let e = [E3[q].0, E3[q].1, E3[q].2];
            let source = core::array::from_fn(|a| (cell[a] as isize - e[a] as isize) as usize);
            post[index(source)][q]
        }
    })
}

#[test]
fn profiled_wall_seams_match_prescribed_force_reference_and_close_transport() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xACC30F);
    let acceleration = [-0.004, 0.003, 0.002];
    for face in Face3::ALL {
        let initial = seamed_grid(face, &gate);
        let reference = prescribed_force_reference(&initial, |_| acceleration, &gate);
        let post: Vec<_> = cells()
            .map(|[x, y, z]| {
                collide_cell3(
                    initial.populations(x, y, z),
                    model(face),
                    initial.force_density(x, y, z),
                )
                .unwrap()
            })
            .collect();
        let links = initial.boundary_links();
        let mut serial = initial.clone();
        let mut pooled = initial.clone();
        let walls = serial
            .step_with_endpoint_body_acceleration_and_wall_momentum(
                acceleration,
                &links,
                [0.2, -0.4, 1.1],
                &gate,
            )
            .unwrap();
        let pooled_walls = pooled
            .step_pooled_with_endpoint_body_acceleration_and_wall_momentum(
                &pool,
                &gate,
                acceleration,
                &links,
                [0.2, -0.4, 1.1],
            )
            .unwrap();
        assert_eq!(snapshot(&serial), snapshot(&pooled));
        assert_eq!(walls, pooled_walls);
        let mut reservoir_mass = 0.0;
        let mut reservoir_impulse = [0.0; 3];
        let axis = face as usize / 2;
        for cell @ [x, y, z] in cells() {
            let actual = serial.populations(x, y, z);
            let expected = reference.populations(x, y, z);
            let (rho, momentum) = moments(&serial, cell);
            for q in 0..Q3 {
                close(actual[q], expected[q], 2e-13);
            }
            for a in 0..3 {
                assert_eq!(serial.force_density(x, y, z)[a], rho * acceleration[a]);
                close(
                    serial.force_density(x, y, z)[a],
                    reference.force_density(x, y, z)[a],
                    2e-14,
                );
            }
            let open = initial.open_link_mask(x, y, z);
            if open == 0 {
                continue;
            }
            let streamed = streamed_before_reconstruction(&initial, &post, cell);
            let wall = initial.link_mask(x, y, z);
            for q in 0..Q3 {
                if wall != 0 && open & (1 << q) == 0 {
                    assert_eq!(
                        actual[q].to_bits(),
                        streamed[q].to_bits(),
                        "face {face:?}, cell {cell:?}, q {q}"
                    );
                }
                let delta = actual[q] - streamed[q];
                let direction = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
                reservoir_mass += delta;
                for a in 0..3 {
                    // The provisional open population is reflected outgoing
                    // data. Count that link's transport, then the actual native
                    // reservoir reconstruction; never call this a closed box.
                    let reflected = if open & (1 << q) != 0 {
                        2.0 * post[index(cell)][OPP3[q]]
                    } else {
                        0.0
                    };
                    reservoir_impulse[a] += (reflected + delta) * direction[a];
                }
            }
            if wall == 0 {
                let mut interior = cell;
                interior[axis] = if cell[axis] == 0 { 1 } else { N - 2 };
                let (neighbor_rho, neighbor_momentum) = moments(&serial, interior);
                let inlet = cell[axis]
                    == if (face as usize).is_multiple_of(2) {
                        0
                    } else {
                        N - 1
                    };
                let (target_rho, target_velocity) = if inlet {
                    (neighbor_rho, velocity_trace(cell, axis))
                } else {
                    (
                        pressure_trace(cell, axis),
                        neighbor_momentum.map(|p| p / neighbor_rho),
                    )
                };
                close(rho, target_rho, 6e-15);
                for a in 0..3 {
                    close(momentum[a] / rho, target_velocity[a], 2e-15);
                }
            }
        }
        let before = totals(&initial);
        let after = totals(&serial);
        assert!(reservoir_mass.abs() > 1e-4);
        close(
            after.0 - before.0,
            walls.fluid_mass_change + reservoir_mass,
            3e-10,
        );
        for a in 0..3 {
            close(
                after.1[a] - before.1[a],
                0.5 * (before.2[a] + after.2[a])
                    + walls.fluid_population_impulse[a]
                    + reservoir_impulse[a],
                4e-11,
            );
        }
    }
}

#[test]
fn spatial_acceleration_uses_each_donor_and_reservoir_cells_own_density() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xACC50F);
    let acceleration = |[x, y, z]: [usize; 3]| {
        [
            0.001 + 0.0002 * x as f64,
            -0.0015 + 0.0003 * y as f64,
            0.0007 - 0.0001 * z as f64,
        ]
    };
    for face in Face3::ALL {
        let initial = seamed_grid(face, &gate);
        let reference = prescribed_force_reference(&initial, acceleration, &gate);
        let mut serial = initial.clone();
        let mut pooled = initial.clone();
        serial
            .step_with_endpoint_body_acceleration_field(acceleration, &gate)
            .unwrap();
        pooled
            .step_pooled_with_endpoint_body_acceleration_field(&pool, &gate, acceleration)
            .unwrap();
        assert_eq!(snapshot(&serial), snapshot(&pooled));
        let axis = face as usize / 2;
        for cell @ [x, y, z] in cells() {
            let actual = serial.populations(x, y, z);
            let expected = reference.populations(x, y, z);
            let (rho, momentum) = moments(&serial, cell);
            for q in 0..Q3 {
                close(actual[q], expected[q], 2e-13);
            }
            assert_eq!(serial.force_density(x, y, z), acceleration(cell).map(|a| rho * a));
            if initial.open_link_mask(x, y, z) == 0 || initial.link_mask(x, y, z) != 0 {
                continue;
            }
            let mut interior = cell;
            interior[axis] = if cell[axis] == 0 { 1 } else { N - 2 };
            let (neighbor_rho, neighbor_momentum) = moments(&serial, interior);
            let inlet = cell[axis]
                == if (face as usize).is_multiple_of(2) {
                    0
                } else {
                    N - 1
                };
            let (target_rho, target_velocity) = if inlet {
                (neighbor_rho, velocity_trace(cell, axis))
            } else {
                (
                    pressure_trace(cell, axis),
                    neighbor_momentum.map(|p| p / neighbor_rho),
                )
            };
            close(rho, target_rho, 6e-15);
            for a in 0..3 {
                close(momentum[a] / rho, target_velocity[a], 2e-15);
            }
        }
    }
}

#[test]
fn failed_open_source_and_initial_bind_preserve_accepted_state_and_allow_recovery() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xACC40F);
    for pooled in [false, true] {
        let mut faces = [FaceBoundary3::stationary_wall(); 6];
        faces[0] = FaceBoundary3::Velocity { velocity: [0.0; 3] };
        faces[1] = FaceBoundary3::Velocity { velocity: [0.0; 3] };
        let mut grid = BoundaryGrid3::new(N, N, N, 0.9, [0.0; 3], BoundarySpec3::new(faces));
        grid.initialize_uniform_equilibrium(2.3, [0.0; 3], &gate)
            .unwrap();
        let before = snapshot(&grid);
        let mut control = grid.clone();
        assert!(
            grid.rebind_body_acceleration([f64::MAX, 0.0, 0.0], &gate)
                .is_err()
        );
        assert_eq!(snapshot(&grid), before);
        assert!(
            RigidFrame3::default()
                .prepare_step(&mut grid, [0.001, 0.0, 0.0], &gate)
                .is_err()
        );
        assert_eq!(snapshot(&grid), before);
        let advance = |grid: &mut BoundaryGrid3, acceleration, gate: &CancelGate| {
            if pooled {
                grid.step_pooled_with_endpoint_body_acceleration(&pool, gate, acceleration)
            } else {
                grid.step_with_endpoint_body_acceleration(acceleration, gate)
            }
        };
        // The min-face closure has already been reconstructed in scratch.
        // At the max-face wall rim, the positive mass response exceeds one;
        // the first-use endpoint field and all proposals must remain private.
        assert_eq!(
            advance(&mut grid, [8.0, 0.0, 0.0], &gate),
            Err(BoundaryStepError3::Unphysical {
                cell: [N - 1, 0, 0]
            }),
        );
        assert_eq!(snapshot(&grid), before);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            grid.rebind_body_acceleration([0.001, 0.0, 0.0], &cancelled)
                .is_err()
        );
        assert_eq!(
            advance(&mut grid, [0.001, 0.0, 0.0], &cancelled),
            Err(BoundaryStepError3::Cancelled)
        );
        assert_eq!(snapshot(&grid), before);
        // A failed initial step must not make topology immutable.
        grid.voxelize_sdf(|_| 1.0);
        grid.rebind_body_acceleration([0.0001, 0.0, 0.0], &gate)
            .unwrap();
        control
            .rebind_body_acceleration([0.0001, 0.0, 0.0], &gate)
            .unwrap();
        advance(&mut grid, [0.0002, 0.0, 0.0], &gate).unwrap();
        advance(&mut control, [0.0002, 0.0, 0.0], &gate).unwrap();
        assert_eq!(snapshot(&grid), snapshot(&control));
    }
}
