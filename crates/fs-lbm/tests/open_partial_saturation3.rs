//! Actual reservoir transport through native single, assembly and dry bodies.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryForcingError3, BoundaryGrid3,
    BoundarySpec3, CollisionModel3, DryBodyCouplingConfig3, E3, Face3, FaceBoundary3, InletFlux3,
    MovingPartialSaturationBody3, PartialSaturationBody3, PartialSaturationError3, Q3, equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const N: usize = 8;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; difference={:.3e}",
        actual - expected
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn models() -> [CollisionModel3; 4] {
    [
        CollisionModel3::Bgk { tau: 0.83 },
        CollisionModel3::Trt {
            tau: 0.83,
            magic: 0.1875,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.0 / 0.83,
            higher_order_rate: 0.93,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.0 / 0.83,
            third_order_rate: 0.97,
            fourth_order_rate: 1.14,
        },
    ]
}

fn boundaries(periodic: bool) -> BoundarySpec3 {
    let wall = FaceBoundary3::stationary_wall();
    let lateral = if periodic {
        FaceBoundary3::Periodic
    } else {
        wall
    };
    BoundarySpec3::new([
        lateral,
        lateral,
        wall,
        wall,
        FaceBoundary3::Velocity {
            velocity: [0.001, -0.002, 0.012],
        },
        FaceBoundary3::Pressure { density: 1.002 },
    ])
}

fn prepared(nx: usize, model: CollisionModel3, boundary: BoundarySpec3) -> BoundaryGrid3 {
    let mut grid = BoundaryGrid3::with_collision_model(nx, N, N, model, [0.0; 3], boundary);
    for z in 0..N {
        for y in 0..N {
            for x in 0..nx {
                let rho = 1.0 + 0.0002 * (x + 2 * y + z) as f64;
                let velocity = [
                    0.0003 * (y as f64 - 3.5),
                    -0.0002 * (x as f64 - 3.5),
                    0.004 + 0.0001 * z as f64,
                ];
                grid.set_populations(x, y, z, &equilibrium3(rho, velocity));
            }
        }
    }
    grid.rebind_force_field(|[x, y, z]| {
        [
            1e-6 * (x as f64 - 3.0),
            -2e-6 * (z as f64 - 2.0),
            1e-6 * (1 + y) as f64,
        ]
    });
    grid
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut values = Vec::with_capacity(nx * ny * nz);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                values.push(grid.populations(x, y, z));
            }
        }
    }
    values
}

// Independent whole-domain population differences. Held source-frame shifts
// cancel here; no response matrix or boundary decomposition enters the sum.
fn fluid_change(before: &[[f64; Q3]], grid: &BoundaryGrid3, origin: [f64; 3]) -> (f64, [f64; 6]) {
    let [nx, ny, nz] = grid.dimensions();
    let mut mass = 0.0;
    let mut change = [0.0; 6];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let old = before[(z * ny + y) * nx + x];
                let after = grid.populations(x, y, z);
                let mut p = [0.0; 3];
                for q in 0..Q3 {
                    let difference = after[q] - old[q];
                    mass += difference;
                    let c = [E3[q].0, E3[q].1, E3[q].2];
                    for i in 0..3 {
                        p[i] += f64::from(c[i]) * difference;
                    }
                }
                let angular = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    p,
                );
                for i in 0..3 {
                    change[i] += p[i];
                    change[i + 3] += angular[i];
                }
            }
        }
    }
    (mass, change)
}

fn body(center: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass: 0.8,
        inertia_world: [[1.3, 0.12, -0.08], [0.12, 1.7, 0.06], [-0.08, 0.06, 2.0]],
        velocity: [-0.004, 0.003, 0.002],
        angular_velocity: [0.002, -0.001, 0.003],
    }
}

fn body_change(
    before: PartialSaturationBody3,
    after: PartialSaturationBody3,
    origin: [f64; 3],
) -> [f64; 6] {
    let p = core::array::from_fn(|i| before.mass * (after.velocity[i] - before.velocity[i]));
    let dw: [f64; 3] =
        core::array::from_fn(|i| after.angular_velocity[i] - before.angular_velocity[i]);
    let spin = before
        .inertia_world
        .map(|row| row.into_iter().zip(dw).map(|(a, b)| a * b).sum::<f64>());
    let angular = cross(core::array::from_fn(|i| before.center[i] - origin[i]), p);
    core::array::from_fn(|i| {
        if i < 3 {
            p[i]
        } else {
            spin[i - 3] + angular[i - 3]
        }
    })
}

#[test]
fn empty_overlay_keeps_native_open_step_bits_for_every_collision_model() {
    for model in models() {
        for periodic in [false, true] {
            let mut native = prepared(N, model, boundaries(periodic));
            let mut coupled = prepared(N, model, boundaries(periodic));
            let mut rigid = body([4.0; 3]);
            let fractions = vec![0.0; N * N * N];
            for _ in 0..3 {
                let before = snapshot(&coupled);
                let origin = [-3.0, 7.0, 2.5];
                native.step();
                let receipt = coupled
                    .step_with_partial_saturation_body(
                        &fractions,
                        &mut rigid,
                        [0.0; 6],
                        origin,
                        &CancelGate::new(),
                    )
                    .unwrap();
                assert_eq!(snapshot(&native), snapshot(&coupled));
                let (mass, change) = fluid_change(&before, &coupled, origin);
                close(mass, receipt.open.mass_change, 2e-12);
                for i in 0..3 {
                    close(
                        change[i],
                        receipt.open.impulse[i]
                            + receipt.walls.fluid_population_impulse[i]
                            + receipt.forcing.impulse[i],
                        4e-12,
                    );
                    close(
                        change[i + 3],
                        receipt.open.angular_impulse[i]
                            + receipt.walls.fluid_population_angular_impulse[i]
                            + receipt.periodic.angular_momentum_change[i]
                            + receipt.forcing.angular_impulse[i],
                        3e-11,
                    );
                }
                assert_eq!(receipt.open.reconstructed_cells, 2 * N * N);
                // Four diagonal directions plus one normal; corners crossing
                // walls belong to walls, including the mixed periodic case.
                assert_eq!(
                    receipt.open.open_links,
                    2 * (5 * N * N - if periodic { 2 * N } else { 4 * N })
                );
            }
        }
    }
}

#[test]
fn simultaneous_bodies_close_forced_reservoir_wall_and_periodic_momentum() {
    for model in models() {
        let origin = [-2.3, 6.7, -1.1];
        let shift = [8.0, -3.0, 2.0];
        let mut first_receipt = None;
        for shifted in [false, true] {
            let origin = core::array::from_fn(|i| origin[i] + if shifted { shift[i] } else { 0.0 });
            let mut grid = prepared(N, model, boundaries(true));
            let before = snapshot(&grid);
            let mut bodies = [body([3.1, 3.4, 3.8]), body([4.8, 4.1, 4.0])];
            let old = bodies;
            let mut fractions = [vec![0.0; N * N * N], vec![0.0; N * N * N]];
            for (x, y, z, a, b) in [
                (3, 3, 3, 0.5, 0.2),
                (4, 4, 3, 0.15, 0.6),
                (2, 4, 4, 0.7, 0.0),
                (5, 3, 4, 0.0, 0.8),
            ] {
                fractions[0][(z * N + y) * N + x] = a;
                fractions[1][(z * N + y) * N + x] = b;
            }
            let external = [
                [0.0003, -0.0002, 0.0001, 0.0001, 0.0004, -0.0002],
                [-0.0002, 0.0004, 0.0003, -0.0002, 0.0001, 0.0003],
            ];
            let receipt = grid
                .step_with_partial_saturation_bodies(
                    fractions.each_ref().map(Vec::as_slice),
                    &mut bodies,
                    external,
                    origin,
                    &CancelGate::new(),
                )
                .unwrap();
            let (mass, mut change) = fluid_change(&before, &grid, origin);
            let mut load = [0.0; 6];
            for a in 0..2 {
                let dp = body_change(old[a], bodies[a], origin);
                let moment = cross(
                    core::array::from_fn(|i| old[a].center[i] - origin[i]),
                    external[a][..3].try_into().unwrap(),
                );
                for i in 0..6 {
                    change[i] += dp[i];
                    load[i] += external[a][i] + if i < 3 { 0.0 } else { moment[i - 3] };
                }
            }
            close(mass, receipt.open.mass_change, 2e-12);
            for i in 0..3 {
                close(
                    change[i],
                    load[i]
                        + receipt.open.impulse[i]
                        + receipt.walls.fluid_population_impulse[i]
                        + receipt.forcing.impulse[i],
                    4e-12,
                );
                close(
                    change[i + 3],
                    load[i + 3]
                        + receipt.open.angular_impulse[i]
                        + receipt.walls.fluid_population_angular_impulse[i]
                        + receipt.periodic.angular_momentum_change[i]
                        + receipt.forcing.angular_impulse[i],
                    3e-11,
                );
            }
            if let Some(original) = first_receipt {
                let original: fs_lbm::d3q19::OpenBoundaryMomentumTransport3 = original;
                assert_eq!(original.impulse, receipt.open.impulse);
                assert_eq!(original.mass_change, receipt.open.mass_change);
                let moment_shift = cross(shift, original.impulse);
                for i in 0..3 {
                    close(
                        receipt.open.angular_impulse[i],
                        original.angular_impulse[i] - moment_shift[i],
                        2e-12,
                    );
                }
            } else {
                first_receipt = Some(receipt.open);
            }
        }
    }
}

#[test]
fn occupied_open_faces_and_cancellation_refuse_without_changing_fluid_or_bodies() {
    let mut grid = prepared(N, models()[2], boundaries(false));
    let before = snapshot(&grid);
    let mut bodies = [body([3.0; 3]), body([4.0; 3])];
    let old = bodies;
    let mut fractions = [vec![0.0; N * N * N], vec![0.0; N * N * N]];
    fractions[1][(2 * N) + 3] = 0.2;
    assert!(matches!(
        grid.step_with_partial_saturation_bodies(
            fractions.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Cell {
            cell: [3, 2, 0],
            ..
        })
    ));
    assert_eq!(snapshot(&grid), before);
    assert_eq!(bodies, old);
    let mut one = old[0];
    assert!(matches!(
        grid.step_with_partial_saturation_body(
            &fractions[1],
            &mut one,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Cell { .. })
    ));
    assert_eq!(one, old[0]);
    fractions[1].fill(0.0);
    fractions[0][(3 * N + 3) * N + 3] = 0.7;
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        grid.step_with_partial_saturation_bodies(
            fractions.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(snapshot(&grid), before);
    assert_eq!(bodies, old);
    let mut clean = prepared(N, models()[2], boundaries(false));
    let mut clean_bodies = old;
    let retry = grid
        .step_with_partial_saturation_bodies(
            fractions.each_ref().map(Vec::as_slice),
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    let reference = clean
        .step_with_partial_saturation_bodies(
            fractions.each_ref().map(Vec::as_slice),
            &mut clean_bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(retry, reference);
    assert_eq!(snapshot(&grid), snapshot(&clean));
    assert_eq!(bodies, clean_bodies);
}

#[test]
fn late_open_operator_overflow_keeps_accepted_state_and_scratch_reusable() {
    let wall = FaceBoundary3::stationary_wall();
    let boundary = BoundarySpec3::new([
        wall,
        wall,
        wall,
        wall,
        FaceBoundary3::Velocity {
            velocity: [0.0, 0.0, 0.01],
        },
        FaceBoundary3::Pressure { density: 1e308 },
    ]);
    let mut grid = BoundaryGrid3::new(N, N, N, 0.83, [0.0; 3], boundary);
    let before = snapshot(&grid);
    let mut rigid = body([4.0; 3]);
    let old = rigid;
    let fractions = vec![0.0; N * N * N];
    for _ in 0..2 {
        // Collision and the first (velocity) face succeed. The opposite
        // reservoir cannot produce a finite total exchange; publication fails.
        assert!(
            grid.step_with_partial_saturation_body(
                &fractions,
                &mut rigid,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(snapshot(&grid), before);
        assert_eq!(rigid, old);
    }
    let mut bodies = [old; 2];
    assert!(
        grid.step_with_partial_saturation_bodies(
            [&fractions; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(snapshot(&grid), before);
    assert_eq!(bodies, [old; 2]);
}

#[test]
fn cancellation_after_open_reconstruction_preserves_moving_body_and_replays() {
    struct CancelAtEndpoint<'a> {
        shape: BodyFrameEllipsoid3,
        gate: &'a CancelGate,
        radius_calls: std::cell::Cell<usize>,
    }
    impl BodyFrameShape3 for CancelAtEndpoint<'_> {
        fn bounding_radius(&self) -> f64 {
            let calls = self.radius_calls.get() + 1;
            self.radius_calls.set(calls);
            // Initial admission, initial sampling, midpoint sampling, then
            // endpoint sampling in the post-collision/publication callback.
            // Geometry is unchanged; only the real cancellation gate changes.
            if calls == 4 {
                self.gate.request();
            }
            self.shape.bounding_radius()
        }
        fn contains(&self, point: Vec3) -> bool {
            self.shape.contains(point)
        }
    }
    let mut grid = prepared(N, models()[1], boundaries(true));
    let before = snapshot(&grid);
    let mut body = MovingPartialSaturationBody3 {
        center: [4.0; 3],
        orientation: So3::identity(),
        mass: 0.8,
        inertia: InertiaTensor3::new(Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        })
        .unwrap(),
        velocity: [0.001, -0.002, 0.003],
        angular_momentum_body: Vec3::new(0.001, -0.001, 0.002),
    };
    let initial = body;
    let gate = CancelGate::new();
    let shape = CancelAtEndpoint {
        shape: BodyFrameEllipsoid3::sphere(1.1).unwrap(),
        gate: &gate,
        radius_calls: std::cell::Cell::new(0),
    };
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    let refused = grid.step_with_moving_partial_saturation_body(
        &mut sampler,
        &shape,
        &mut body,
        [0.0; 6],
        [0.0; 3],
        &gate,
    );
    assert!(matches!(refused, Err(PartialSaturationError3::Cancelled)));
    assert_eq!(shape.radius_calls.get(), 4);
    assert_eq!(body, initial);
    assert_eq!(snapshot(&grid), before);
    let mut fresh = prepared(N, models()[1], boundaries(true));
    let mut fresh_sampler = BodyFrameSampler3::new(&fresh, 4).unwrap();
    let mut fresh_body = initial;
    let retry = grid
        .step_with_moving_partial_saturation_body(
            &mut sampler,
            &shape.shape,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    let reference = fresh
        .step_with_moving_partial_saturation_body(
            &mut fresh_sampler,
            &shape.shape,
            &mut fresh_body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(retry, reference);
    assert_eq!(snapshot(&grid), snapshot(&fresh));
    assert_eq!(body, fresh_body);
}

// Direct physical exterior-fluid-plus-dry-solid inventory. The native sampler
// exposes the declared cut-cell first moments; no fictitious interior momentum
// or pressure-reservoir load is assigned to the dry body's mass.
fn dry_inventory<const K: usize>(
    grid: &BoundaryGrid3,
    samplers: &[BodyFrameSampler3; K],
    bodies: [MovingPartialSaturationBody3; K],
    origin: [f64; 3],
) -> [f64; 6] {
    let [nx, ny, nz] = grid.dimensions();
    let mut first = vec![[0.0; 3]; nx * ny * nz];
    for sampler in samplers {
        for &(index, moment) in sampler.cut_cell_first_moments() {
            for i in 0..3 {
                first[index][i] += moment[i];
            }
        }
    }
    let mut result = [0.0; 6];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let index = (z * ny + y) * nx + x;
                let exterior = 1.0 - samplers.iter().map(|s| s.fractions()[index]).sum::<f64>();
                let populations = grid.populations(x, y, z);
                let force = grid.force_density(x, y, z);
                let p: [f64; 3] = core::array::from_fn(|i| {
                    0.5 * force[i]
                        + (0..Q3)
                            .map(|q| {
                                let c = [E3[q].0, E3[q].1, E3[q].2];
                                populations[q] * f64::from(c[i])
                            })
                            .sum::<f64>()
                });
                let angular = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    p,
                );
                let excluded = cross(first[index], p);
                for i in 0..3 {
                    result[i] += exterior * p[i];
                    result[i + 3] += exterior * angular[i] - excluded[i];
                }
            }
        }
    }
    for body in bodies {
        let p = body.velocity.map(|v| body.mass * v);
        let spin = body.orientation.rotate(body.angular_momentum_body).unwrap();
        let spin = [spin.x, spin.y, spin.z];
        let angular = cross(core::array::from_fn(|i| body.center[i] - origin[i]), p);
        for i in 0..3 {
            result[i] += p[i];
            result[i + 3] += angular[i] + spin[i];
        }
    }
    result
}

fn dry_case<const K: usize>(single: bool) {
    let mut grid = prepared(12, models()[0], boundaries(true));
    let gate = CancelGate::new();
    let shape = BodyFrameEllipsoid3::sphere(1.1).unwrap();
    let mut samplers: [_; K] = core::array::from_fn(|_| BodyFrameSampler3::new(&grid, 4).unwrap());
    let mut bodies = core::array::from_fn(|a| {
        let mass = (0.6 + a as f64) * 4.0 * core::f64::consts::PI * 1.1_f64.powi(3) / 3.0;
        MovingPartialSaturationBody3 {
            center: [3.0 + 3.0 * a as f64, 4.0, 4.0],
            orientation: So3::identity(),
            mass,
            inertia: InertiaTensor3::new(Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
                    .map(|v| 0.4 * mass * 1.1_f64.powi(2) * v),
            })
            .unwrap(),
            velocity: [0.0; 3],
            angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
        }
    });
    let origin = [-2.0, 6.0, -3.0];
    for a in 0..K {
        samplers[a]
            .sample(bodies[a].center, bodies[a].orientation, &shape, &gate)
            .unwrap();
    }
    let before = dry_inventory(&grid, &samplers, bodies, origin);
    let (open, walls, periodic, forcing, residual) = if single {
        assert_eq!(K, 1);
        let receipt = grid
            .step_with_dry_partial_saturation_body(
                &mut samplers[0],
                &shape,
                &mut bodies[0],
                [0.0; 6],
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let native = receipt.native.coupling;
        (
            native.open,
            native.walls,
            native.periodic,
            receipt.represented_fluid_forcing_impulse,
            receipt.combined_momentum_residual,
        )
    } else {
        let receipt = grid
            .step_with_dry_partial_saturation_bodies(
                &mut samplers,
                [&shape; K],
                &mut bodies,
                [[0.0; 6]; K],
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let native = receipt.native.coupling;
        for response in receipt.native.bodies {
            assert_eq!(response.coupling.open, Default::default());
        }
        (
            native.open,
            native.walls,
            native.periodic,
            receipt.represented_fluid_forcing_impulse,
            receipt.combined_momentum_residual,
        )
    };
    for a in 0..K {
        samplers[a]
            .sample(bodies[a].center, bodies[a].orientation, &shape, &gate)
            .unwrap();
    }
    let after = dry_inventory(&grid, &samplers, bodies, origin);
    assert!(open.impulse.iter().any(|v| v.abs() > 1e-3));
    for i in 0..3 {
        close(
            after[i] - before[i],
            open.impulse[i] + walls.fluid_population_impulse[i] + forcing[i],
            3e-10,
        );
        close(
            after[i + 3] - before[i + 3],
            open.angular_impulse[i]
                + walls.fluid_population_angular_impulse[i]
                + periodic.angular_momentum_change[i]
                + forcing[i + 3],
            2e-9,
        );
    }
    for value in residual {
        close(value, 0.0, 3e-10);
    }
}

#[test]
fn physical_dry_single_singleton_and_pair_include_reservoir_momentum_once() {
    dry_case::<1>(true);
    dry_case::<1>(false);
    dry_case::<2>(false);
}

#[test]
fn live_targets_reuse_warmed_forced_solver_and_measure_actual_reservoir_transport() {
    let gate = CancelGate::new_clock_free();
    let origin = [-2.3, 6.7, -1.1];
    let mut fractions = vec![0.0; N * N * N];
    fractions[(3 * N + 3) * N + 3] = 0.45;
    for model in models() {
        let mut actual = prepared(N, model, boundaries(true));
        let mut rigid = body([4.0; 3]);
        actual
            .step_with_partial_saturation_body(
                &fractions, &mut rigid, [0.0; 6], origin, &gate,
            )
            .unwrap();
        let forces: Vec<_> = (0..N * N * N)
            .map(|i| actual.force_density(i % N, (i / N) % N, i / (N * N)))
            .collect();
        for (velocity, density) in [
            ([0.002, -0.001, 0.006], 1.001),
            ([-0.001, 0.001, -0.004], 0.999),
            ([0.001, -0.002, 0.012], 1.002),
        ] {
            let before = snapshot(&actual);
            let old_body = rigid;
            let mut faces = Face3::ALL.map(|face| actual.boundaries().face(face));
            faces[4] = FaceBoundary3::Velocity { velocity };
            faces[5] = FaceBoundary3::Pressure { density };
            let targets = BoundarySpec3::new(faces);
            actual.rebind_open_boundaries(targets, &gate).unwrap();
            assert_eq!(actual.boundaries(), targets);
            assert_eq!(snapshot(&actual), before);
            for (i, &force) in forces.iter().enumerate() {
                assert_eq!(actual.force_density(i % N, (i / N) % N, i / (N * N)), force);
            }
            // Independent construction with the new controls and the same
            // accepted populations/force frame, rather than another live rebind.
            let mut reference = prepared(N, model, targets);
            for (i, values) in before.iter().enumerate() {
                reference.set_populations(i % N, (i / N) % N, i / (N * N), values);
            }
            let mut reference_body = old_body;
            let expected = reference
                .step_with_partial_saturation_body(
                    &fractions, &mut reference_body, [0.0; 6], origin, &gate,
                )
                .unwrap();
            let receipt = actual
                .step_with_partial_saturation_body(
                    &fractions, &mut rigid, [0.0; 6], origin, &gate,
                )
                .unwrap();
            assert_eq!(receipt, expected);
            assert_eq!(rigid, reference_body);
            assert_eq!(snapshot(&actual), snapshot(&reference));
            let (mass, fluid) = fluid_change(&before, &actual, origin);
            let body = body_change(old_body, rigid, origin);
            close(mass, receipt.open.mass_change, 2e-12);
            for i in 0..3 {
                close(
                    fluid[i] + body[i],
                    receipt.open.impulse[i]
                        + receipt.walls.fluid_population_impulse[i]
                        + receipt.forcing.impulse[i],
                    4e-12,
                );
                close(
                    fluid[i + 3] + body[i + 3],
                    receipt.open.angular_impulse[i]
                        + receipt.walls.fluid_population_angular_impulse[i]
                        + receipt.periodic.angular_momentum_change[i]
                        + receipt.forcing.angular_impulse[i],
                    3e-11,
                );
            }
            assert_eq!(receipt.open.reconstructed_cells, 2 * N * N);
            assert!(receipt.open.impulse.iter().any(|value| value.abs() > 1e-6));
        }
    }
}

#[test]
fn live_target_refusals_preserve_profiles_populations_and_previous_controls() {
    let gate = CancelGate::new_clock_free();
    let original = boundaries(true);
    let mut grid = prepared(N, models()[0], original);
    grid.set_velocity_profile(Face3::ZMin, InletFlux3::Unconstrained, &gate, |_| {
        Ok([0.001, -0.002, 0.012])
    })
    .unwrap();
    let accepted = snapshot(&grid);
    let profile = grid.velocity_profile(Face3::ZMin).unwrap().to_vec();
    let original_faces = Face3::ALL.map(|face| original.face(face));
    let mut proposals = Vec::new();
    let mut changed = original_faces;
    changed[0] = FaceBoundary3::stationary_wall();
    changed[1] = FaceBoundary3::stationary_wall();
    proposals.push(BoundarySpec3::new(changed));
    changed = original_faces;
    changed[4] = FaceBoundary3::Pressure { density: 1.0 };
    proposals.push(BoundarySpec3::new(changed));
    changed = original_faces;
    changed[2] = FaceBoundary3::Wall { velocity: [0.001, 0.0, 0.0] };
    proposals.push(BoundarySpec3::new(changed));
    for velocity in [[f64::NAN, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.0, 0.02]] {
        changed = original_faces;
        changed[4] = FaceBoundary3::Velocity { velocity };
        changed[5] = FaceBoundary3::Pressure { density: 1.001 };
        proposals.push(BoundarySpec3::new(changed));
    }
    for density in [0.0, f64::NAN, f64::INFINITY] {
        changed = original_faces;
        changed[5] = FaceBoundary3::Pressure { density };
        proposals.push(BoundarySpec3::new(changed));
    }
    for targets in proposals {
        assert!(matches!(
            grid.rebind_open_boundaries(targets, &gate),
            Err(BoundaryForcingError3::Invalid(_))
        ));
        assert_eq!(grid.boundaries(), original);
        assert_eq!(snapshot(&grid), accepted);
        assert_eq!(grid.velocity_profile(Face3::ZMin).unwrap(), profile.as_slice());
    }
    changed = original_faces;
    changed[5] = FaceBoundary3::Pressure { density: 1.001 };
    let pressure_only = BoundarySpec3::new(changed);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        grid.rebind_open_boundaries(pressure_only, &cancelled),
        Err(BoundaryForcingError3::Cancelled)
    );
    assert_eq!(grid.boundaries(), original);
    assert_eq!(snapshot(&grid), accepted);
    grid.rebind_open_boundaries(pressure_only, &gate).unwrap();
    assert_eq!(grid.velocity_profile(Face3::ZMin).unwrap(), profile.as_slice());
    grid.clear_velocity_profile(Face3::ZMin);
    changed[4] = FaceBoundary3::Velocity { velocity: [0.0, 0.0, 0.02] };
    grid.rebind_open_boundaries(BoundarySpec3::new(changed), &gate).unwrap();
    assert_eq!(snapshot(&grid), accepted);
}
