//! G0/G3: native 3D fluid/body momentum, work, covariance and refusal.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, CollisionModel3, E3, FaceBoundary3, PartialSaturationBody3,
    PartialSaturationError3, Q3, equilibrium3,
};

const N: usize = 8;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; difference={:.4e}",
        a - b
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn body() -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [3.2, 3.7, 3.4],
        mass: 0.7,
        inertia_world: [[1.3, 0.12, -0.08], [0.12, 1.7, 0.06], [-0.08, 0.06, 2.0]],
        velocity: [-0.014, 0.008, 0.006],
        angular_velocity: [0.004, -0.003, 0.005],
    }
}

fn grid(boundaries: BoundarySpec3) -> BoundaryGrid3 {
    BoundaryGrid3::new(N, N, N, 0.83, [0.0; 3], boundaries)
}

fn index(cell: [usize; 3]) -> usize {
    (cell[2] * N + cell[1]) * N + cell[0]
}

fn fractions() -> Vec<f64> {
    let mut fractions = vec![0.0; N * N * N];
    for (cell, epsilon) in [
        ([4, 3, 3], 0.7),
        ([3, 5, 3], 0.45),
        ([2, 3, 4], 1.0),
        ([3, 2, 2], 0.2),
    ] {
        fractions[index(cell)] = epsilon;
    }
    fractions
}

fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut values = Vec::with_capacity(N * N * N);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                values.push(if grid.is_solid(x, y, z) {
                    [0.0; Q3]
                } else {
                    grid.populations(x, y, z)
                });
            }
        }
    }
    values
}

// Direct population difference, independent of the response matrix and receipt.
fn fluid_change(
    before: &[[f64; Q3]],
    after: &BoundaryGrid3,
    origin: [f64; 3],
) -> (f64, [f64; 3], [f64; 3]) {
    let mut mass = 0.0;
    let mut p = [0.0; 3];
    let mut h = [0.0; 3];
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let f = after.populations(x, y, z);
                let old = before[index([x, y, z])];
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    let delta = f[q] - old[q];
                    mass += delta;
                    local[0] += f64::from(E3[q].0) * delta;
                    local[1] += f64::from(E3[q].1) * delta;
                    local[2] += f64::from(E3[q].2) * delta;
                }
                let r = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                let angular = cross(r, local);
                for d in 0..3 {
                    p[d] += local[d];
                    h[d] += angular[d];
                }
            }
        }
    }
    (mass, p, h)
}

fn body_change(
    before: PartialSaturationBody3,
    after: PartialSaturationBody3,
    origin: [f64; 3],
) -> ([f64; 3], [f64; 3], f64) {
    let p = core::array::from_fn(|d| before.mass * (after.velocity[d] - before.velocity[d]));
    let dw = core::array::from_fn(|d| after.angular_velocity[d] - before.angular_velocity[d]);
    let spin = before.inertia_world.map(|row| dot(row, dw));
    let arm = core::array::from_fn(|d| before.center[d] - origin[d]);
    let orbital = cross(arm, p);
    let h = core::array::from_fn(|d| spin[d] + orbital[d]);
    let mean_v = core::array::from_fn(|d| 0.5 * (before.velocity[d] + after.velocity[d]));
    let mean_w =
        core::array::from_fn(|d| 0.5 * (before.angular_velocity[d] + after.angular_velocity[d]));
    (p, h, dot(p, mean_v) + dot(spin, mean_w))
}

#[test]
fn off_center_six_dof_body_closes_actual_streamed_momentum_and_work() {
    let mut grid = grid(BoundarySpec3::periodic());
    grid.set_populations(4, 3, 3, &equilibrium3(1.07, [0.024, -0.017, 0.011]));
    grid.set_populations(2, 3, 4, &equilibrium3(0.96, [-0.012, 0.022, -0.008]));
    let populations = snapshot(&grid);
    let mut body = body();
    let before = body;
    let external = [0.0005, -0.0003, 0.0002, -0.0004, 0.0006, 0.0003];
    let origin = [-4.2, 7.1, 2.3];
    let receipt = grid
        .step_with_partial_saturation_body(
            &fractions(),
            &mut body,
            external,
            origin,
            &CancelGate::new(),
        )
        .unwrap();
    let (mass, fluid_p, fluid_h) = fluid_change(&populations, &grid, origin);
    let (body_p, body_h, energy) = body_change(before, body, origin);
    let arm = core::array::from_fn(|d| before.center[d] - origin[d]);
    let orbital = cross(arm, external[..3].try_into().unwrap());
    close(mass, 0.0, 3e-13);
    for d in 0..3 {
        close(fluid_p[d] + body_p[d], external[d], 2e-13);
        close(fluid_h[d] + body_h[d], external[d + 3] + orbital[d], 2e-12);
        close(
            receipt.hydrodynamic.fluid_population_impulse[d],
            -receipt.hydrodynamic.body_impulse[d],
            0.0,
        );
        close(
            receipt.hydrodynamic.fluid_population_angular_impulse[d],
            -receipt.hydrodynamic.body_angular_impulse[d],
            0.0,
        );
    }
    assert_eq!(receipt.hydrodynamic.coupled_cells, 4);
    assert_eq!(receipt.walls.measured_links, 0);
    for defect in receipt.momentum_residual {
        close(defect, 0.0, 2e-14);
    }
    close(energy, receipt.kinetic_energy_change, 2e-17);
    close(
        energy,
        receipt.hydrodynamic.body_work + receipt.external_work,
        2e-16,
    );
    assert_ne!(body.angular_velocity, before.angular_velocity);
    assert_eq!(body.center, before.center);
    assert_eq!(body.inertia_world, before.inertia_world);
}

#[test]
fn six_resolved_modes_keep_extremely_light_bodies_bounded() {
    let mut grid = grid(BoundarySpec3::periodic());
    let mut fractions = vec![0.0; N * N * N];
    for cell in [
        [2, 3, 3],
        [4, 3, 3],
        [3, 2, 3],
        [3, 4, 3],
        [3, 3, 2],
        [3, 3, 4],
    ] {
        fractions[index(cell)] = 0.5;
    }
    let initial_v = [0.012, -0.015, 0.007];
    let initial_w = [0.003, -0.002, 0.004];
    let mut body = PartialSaturationBody3 {
        center: [3.5; 3],
        mass: 1e-18,
        inertia_world: [[1e-18, 0.0, 0.0], [0.0, 2e-18, 0.0], [0.0, 0.0, 3e-18]],
        velocity: initial_v,
        angular_velocity: initial_w,
    };
    let receipt = grid
        .step_with_partial_saturation_body(
            &fractions,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    // Native midpoint damping tends to a bounded sign flip as M/K -> 0,
    // unlike an explicit force divided by the tiny mass.
    for d in 0..3 {
        close(body.velocity[d], -initial_v[d], 2e-15);
        close(body.angular_velocity[d], -initial_w[d], 2e-15);
    }
    for defect in receipt.momentum_residual {
        close(defect, 0.0, 5e-16);
    }
}

#[test]
fn common_translation_and_zero_fraction_recover_native_fluid() {
    let mut uniform = grid(BoundarySpec3::periodic());
    let velocity = [0.025, -0.018, 0.013];
    let equilibrium = equilibrium3(1.02, velocity);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                uniform.set_populations(x, y, z, &equilibrium);
            }
        }
    }
    let mut body = body();
    body.velocity = velocity;
    body.angular_velocity = [0.0; 3];
    for _ in 0..16 {
        let receipt = uniform
            .step_with_partial_saturation_body(
                &fractions(),
                &mut body,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        for load in receipt.hydrodynamic.body_impulse {
            close(load, 0.0, 2e-14);
        }
    }
    for d in 0..3 {
        close(body.velocity[d], velocity[d], 2e-14);
        close(body.angular_velocity[d], 0.0, 2e-14);
    }
    for f in snapshot(&uniform) {
        for q in 0..Q3 {
            close(f[q], equilibrium[q], 2e-14);
        }
    }

    let boundaries = BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]);
    let mut reference = grid(boundaries);
    reference.perturb(712, 0.01);
    let mut coupled = grid(boundaries);
    coupled.perturb(712, 0.01);
    let zeros = vec![0.0; N * N * N];
    let mut rest = body;
    rest.velocity = [0.0; 3];
    rest.angular_velocity = [0.0; 3];
    for _ in 0..8 {
        reference.step();
        let receipt = coupled
            .step_with_partial_saturation_body(
                &zeros,
                &mut rest,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(receipt.hydrodynamic.coupled_cells, 0);
        assert_eq!(snapshot(&reference), snapshot(&coupled));
    }
}

#[test]
fn stationary_wall_exchange_closes_repeated_full_domain_balance() {
    let mut grid = grid(BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]));
    let mut body = body();
    let fractions = fractions();
    let origin = [-2.0, 1.3, 5.7];
    let gate = CancelGate::new();
    let mut wall_activity = 0.0;
    for step in 0..20 {
        let old = snapshot(&grid);
        let before = body;
        let external = [1e-4, 0.0, -2e-4, 0.0, 3e-5, -4e-5];
        let receipt = grid
            .step_with_partial_saturation_body(&fractions, &mut body, external, origin, &gate)
            .unwrap();
        let (mass, p, h) = fluid_change(&old, &grid, origin);
        let (bp, bh, _) = body_change(before, body, origin);
        let arm = core::array::from_fn(|d| before.center[d] - origin[d]);
        let ext_h = cross(arm, external[..3].try_into().unwrap());
        close(mass, 0.0, 4e-13);
        for d in 0..3 {
            close(
                p[d] + bp[d] + receipt.walls.wall_impulse[d],
                external[d],
                1e-11,
            );
            close(
                h[d] + bh[d] + receipt.walls.wall_angular_impulse[d],
                external[d + 3] + ext_h[d],
                2e-10,
            );
            if step > 5 {
                wall_activity += receipt.walls.wall_impulse[d].abs();
            }
        }
        close(
            receipt.kinetic_energy_change,
            receipt.hydrodynamic.body_work + receipt.external_work,
            3e-16,
        );
    }
    assert!(wall_activity > 1e-4, "wall accounting must be exercised");
}

fn cycle(v: [f64; 3]) -> [f64; 3] {
    [v[2], v[0], v[1]]
}
fn cycle_cell(v: [usize; 3]) -> [usize; 3] {
    [v[2], v[0], v[1]]
}

#[test]
fn proper_axis_rotation_and_moment_origin_changes_are_covariant() {
    let mut first = grid(BoundarySpec3::periodic());
    let mut second = grid(BoundarySpec3::periodic());
    let cell = [4, 3, 3];
    let velocity = [0.031, -0.017, 0.009];
    first.set_populations(cell[0], cell[1], cell[2], &equilibrium3(1.03, velocity));
    let c = cycle_cell(cell);
    second.set_populations(c[0], c[1], c[2], &equilibrium3(1.03, cycle(velocity)));
    let a = fractions();
    let mut b = vec![0.0; N * N * N];
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                b[index(cycle_cell([x, y, z]))] = a[index([x, y, z])];
            }
        }
    }
    let mut body_a = body();
    let mut body_b = body_a;
    body_b.center = cycle(body_a.center);
    body_b.velocity = cycle(body_a.velocity);
    body_b.angular_velocity = cycle(body_a.angular_velocity);
    let permutation = [2, 0, 1];
    body_b.inertia_world = core::array::from_fn(|i| {
        core::array::from_fn(|j| body_a.inertia_world[permutation[i]][permutation[j]])
    });
    let origin = [-3.2, 1.1, 6.3];
    let shift = [4.0, -2.0, 1.0];
    let origin_b = cycle(core::array::from_fn(|d| origin[d] + shift[d]));
    let receipt_a = first
        .step_with_partial_saturation_body(&a, &mut body_a, [0.0; 6], origin, &CancelGate::new())
        .unwrap();
    let receipt_b = second
        .step_with_partial_saturation_body(&b, &mut body_b, [0.0; 6], origin_b, &CancelGate::new())
        .unwrap();
    let transformed = cross(shift, receipt_a.hydrodynamic.body_impulse);
    let expected_torque = cycle(core::array::from_fn(|d| {
        receipt_a.hydrodynamic.body_angular_impulse[d] - transformed[d]
    }));
    for d in 0..3 {
        close(body_b.velocity[d], cycle(body_a.velocity)[d], 3e-15);
        close(
            body_b.angular_velocity[d],
            cycle(body_a.angular_velocity)[d],
            3e-15,
        );
        close(
            receipt_b.hydrodynamic.body_angular_impulse[d],
            expected_torque[d],
            5e-15,
        );
    }
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let c = cycle_cell([x, y, z]);
                close(
                    second.density(c[0], c[1], c[2]),
                    first.density(x, y, z),
                    3e-15,
                );
                for d in 0..3 {
                    close(
                        second.velocity(c[0], c[1], c[2])[d],
                        cycle(first.velocity(x, y, z))[d],
                        3e-15,
                    );
                }
            }
        }
    }
}

#[test]
fn invalid_requests_and_cancellation_publish_nothing_and_recover() {
    let mut grid = grid(BoundarySpec3::periodic());
    let mut body = body();
    let before = body;
    let populations = snapshot(&grid);
    let good = fractions();
    let mut bad = good.clone();
    bad[N * N * N - 1] = f64::NAN;
    assert!(matches!(
        grid.step_with_partial_saturation_body(
            &bad,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Cell { .. })
    ));
    assert_eq!(body, before);
    assert_eq!(snapshot(&grid), populations);
    // Native collision refusal after earlier cells have written scratch.
    let last = grid.populations(N - 1, N - 1, N - 1);
    grid.set_populations(N - 1, N - 1, N - 1, &[0.0; Q3]);
    let invalid_fluid = snapshot(&grid);
    assert!(matches!(
        grid.step_with_partial_saturation_body(
            &good,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Collision { .. })
    ));
    assert_eq!(body, before);
    assert_eq!(snapshot(&grid), invalid_fluid);
    grid.set_populations(N - 1, N - 1, N - 1, &last);
    body.inertia_world[0][1] += 0.1;
    let bad_body = body;
    assert!(matches!(
        grid.step_with_partial_saturation_body(
            &good,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Invalid { .. })
    ));
    assert_eq!(body, bad_body);
    assert_eq!(snapshot(&grid), populations);
    body = before;
    body.inertia_world[0][0] = -1.0;
    let bad_body = body;
    assert!(
        grid.step_with_partial_saturation_body(
            &good,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(body, bad_body);
    assert_eq!(snapshot(&grid), populations);
    body = before;
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        grid.step_with_partial_saturation_body(&good, &mut body, [0.0; 6], [0.0; 3], &gate),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(body, before);
    assert_eq!(snapshot(&grid), populations);
    // A late endpoint admission failure may write scratch, never populations.
    let mut endpoint = good.clone();
    endpoint[index([7, 7, 7])] = 0.01;
    assert!(
        grid.step_with_partial_saturation_body(
            &endpoint,
            &mut body,
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0; 3],
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(body, before);
    assert_eq!(snapshot(&grid), populations);
    let mut fresh = BoundaryGrid3::new(N, N, N, 0.83, [0.0; 3], BoundarySpec3::periodic());
    let mut fresh_body = before;
    let a = grid
        .step_with_partial_saturation_body(&good, &mut body, [0.0; 6], [0.0; 3], &CancelGate::new())
        .unwrap();
    let b = fresh
        .step_with_partial_saturation_body(
            &good,
            &mut fresh_body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(body, fresh_body);
    assert_eq!(snapshot(&grid), snapshot(&fresh));
}

#[test]
fn unsupported_collision_forcing_open_faces_and_voxel_occupancy_refuse() {
    let mut body = body();
    let good = fractions();
    let mut grids = [
        BoundaryGrid3::with_collision_model(
            N,
            N,
            N,
            CollisionModel3::CentralMoment {
                second_order_rate: 1.2,
                higher_order_rate: 1.5,
            },
            [0.0; 3],
            BoundarySpec3::periodic(),
        ),
        BoundaryGrid3::new(N, N, N, 0.83, [1e-5, 0.0, 0.0], BoundarySpec3::periodic()),
        BoundaryGrid3::new(
            N,
            N,
            N,
            0.83,
            [0.0; 3],
            BoundarySpec3::pressure_duct_z(1.001, 1.0),
        ),
    ];
    for grid in &mut grids {
        let old = snapshot(grid);
        let old_body = body;
        assert!(
            grid.step_with_partial_saturation_body(
                &good,
                &mut body,
                [0.0; 6],
                [0.0; 3],
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(snapshot(grid), old);
        assert_eq!(body, old_body);
    }
    let mut grid = grid(BoundarySpec3::periodic());
    grid.voxelize_sdf(|p| if p == [4.5, 3.5, 3.5] { -1.0 } else { 1.0 });
    let old = snapshot(&grid);
    assert!(matches!(
        grid.step_with_partial_saturation_body(
            &good,
            &mut body,
            [0.0; 6],
            [0.0; 3],
            &CancelGate::new()
        ),
        Err(PartialSaturationError3::Cell { .. })
    ));
    assert_eq!(snapshot(&grid), old);
}
