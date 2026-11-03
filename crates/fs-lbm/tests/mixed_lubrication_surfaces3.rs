//! G1/G3 mixed surface work and momentum; G4 common native publication.
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, CollisionModel3, E3,
    FiniteLubricationPlane3, LubricationPlane3, LubricationSurface3, PartialSaturationBody3, Q3,
    SphereLubrication3,
};

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
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
fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        velocity,
        inertia_world: [
            [0.12, 0.009, -0.004],
            [0.009, 0.17, 0.006],
            [-0.004, 0.006, 0.21],
        ]
        .map(|row| row.map(|x| mass * x)),
        angular_velocity: [0.003, -0.004, 0.002],
    }
}
fn grid(collision: CollisionModel3) -> BoundaryGrid3 {
    BoundaryGrid3::with_collision_model(8, 8, 8, collision, [0.0; 3], BoundarySpec3::periodic())
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.11, 0.1, 0.001, 2)
        .unwrap()
        .with_resistance_mode_limit(10)
        .unwrap()
}
fn sphere() -> PartialSaturationBody3 {
    body([3.5, 4.5, 4.5], 3.0, [-0.008, 0.004, -0.002])
}
fn carrier() -> FiniteLubricationPlane3 {
    FiniteLubricationPlane3::new(
        [2.79, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        body([2.6, 4.1, 4.8], 5.0, [0.002, -0.001, 0.003]),
        0.03,
        0.001,
    )
    .unwrap()
    .with_tangential(0.013, 0.001)
    .unwrap()
}
fn prescribed() -> LubricationPlane3 {
    LubricationPlane3::new(
        [4.21, 4.3, 4.6],
        [-1.0, 0.0, 0.0],
        [0.001, -0.002, 0.001],
        0.03,
        0.001,
    )
    .unwrap()
    .with_tangential(0.013, 0.001)
    .unwrap()
    .with_axial_spin(0.003)
    .unwrap()
}
fn surfaces() -> [LubricationSurface3; 2] {
    [
        LubricationSurface3::Finite(carrier()),
        LubricationSurface3::Prescribed(prescribed()),
    ]
}
fn finite(surface: LubricationSurface3) -> FiniteLubricationPlane3 {
    let LubricationSurface3::Finite(value) = surface else {
        panic!("finite test surface")
    };
    value
}
fn fractions(g: &BoundaryGrid3) -> Vec<f64> {
    let mut sampler = BodyFrameSampler3::new(g, 4).unwrap();
    sampler
        .sample(
            sphere().center,
            So3::identity(),
            &BodyFrameEllipsoid3::sphere(0.7).unwrap(),
            &CancelGate::new(),
        )
        .unwrap();
    sampler.fractions().to_vec()
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..512)
        .map(|i| g.populations(i % 8, (i / 8) % 8, i / 64))
        .collect()
}
fn energy(b: PartialSaturationBody3) -> f64 {
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                b.angular_velocity,
                b.inertia_world.map(|row| dot(row, b.angular_velocity)),
            )
}
fn angular(b: PartialSaturationBody3) -> [f64; 3] {
    let orbital = cross(
        core::array::from_fn(|i| b.center[i] - ORIGIN[i]),
        b.velocity.map(|v| b.mass * v),
    );
    let spin = b.inertia_world.map(|row| dot(row, b.angular_velocity));
    core::array::from_fn(|i| orbital[i] + spin[i])
}
fn external_moment(b: PartialSaturationBody3, j: [f64; 6]) -> [f64; 3] {
    let orbital = cross(
        core::array::from_fn(|i| b.center[i] - ORIGIN[i]),
        core::array::from_fn(|i| j[i]),
    );
    core::array::from_fn(|i| orbital[i] + j[i + 3])
}
fn work(before: PartialSaturationBody3, after: PartialSaturationBody3, j: [f64; 6]) -> f64 {
    (0..3)
        .map(|i| {
            0.5 * j[i] * (before.velocity[i] + after.velocity[i])
                + 0.5 * j[i + 3] * (before.angular_velocity[i] + after.angular_velocity[i])
        })
        .sum()
}

#[test]
fn mixed_surface_recoil_rotation_and_forced_fluid_close_each_shared_tick() {
    for collision in [
        CollisionModel3::Bgk { tau: 0.83 },
        CollisionModel3::Trt {
            tau: 0.83,
            magic: 0.19,
        },
    ] {
        let mut g = grid(collision);
        let force = [2e-7, -1e-7, 0.5e-7];
        g.set_force_field(|cell| {
            if cell.into_iter().all(|i| (2..6).contains(&i)) {
                force
            } else {
                [0.0; 3]
            }
        });
        let fraction = fractions(&g);
        let mut spheres = [sphere()];
        let mut surfaces = surfaces();
        let load = [0.0002, -0.0001, 0.0003, 0.00001, -0.00002, 0.00001];
        let wall_load = [-0.0003, 0.0002, -0.0001, -0.00002, 0.00001, 0.00003];
        for _ in 0..2 {
            let before = snapshot(&g);
            let old_sphere = spheres[0];
            let old_carrier = finite(surfaces[0]).body();
            let report = g
                .step_with_lubricated_partial_saturation_surfaces(
                    [&fraction],
                    [0.7],
                    &mut spheres,
                    [load],
                    model(),
                    &mut surfaces,
                    &[wall_load, [0.0; 6]],
                    ORIGIN,
                    &CancelGate::new(),
                )
                .unwrap();
            let new_carrier = finite(surfaces[0]).body();
            assert_eq!(surfaces[1], LubricationSurface3::Prescribed(prescribed()));
            assert_eq!(
                report.planes[0].finite_body.unwrap().body_after,
                new_carrier
            );
            assert!(report.planes[1].finite_body.is_none());
            assert_eq!(report.coupling.lubrication.active_plane_pairs, 2);
            assert_eq!(report.coupling.lubrication.active_resistance_modes, 10);
            assert!(report.coupling.bodies[0].hydrodynamic.coupled_cells > 0);
            let after = snapshot(&g);
            let mut dp = [0.0; 3];
            let mut dl = [0.0; 3];
            for cell_index in 0..512 {
                let cell = [cell_index % 8, (cell_index / 8) % 8, cell_index / 64];
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    let delta = after[cell_index][q] - before[cell_index][q];
                    for (axis, velocity) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[axis] += f64::from(velocity) * delta;
                    }
                }
                if cell.into_iter().all(|i| (2..6).contains(&i)) {
                    for i in 0..3 {
                        local[i] -= force[i];
                    }
                }
                let torque = cross(
                    core::array::from_fn(|i| cell[i] as f64 + 0.5 - ORIGIN[i]),
                    local,
                );
                for i in 0..3 {
                    dp[i] += local[i];
                    dl[i] += torque[i];
                }
            }
            for (before_body, after_body, external) in [
                (old_sphere, spheres[0], load),
                (old_carrier, new_carrier, wall_load),
            ] {
                let old_l = angular(before_body);
                let new_l = angular(after_body);
                let applied_l = external_moment(before_body, external);
                for i in 0..3 {
                    dp[i] += before_body.mass * (after_body.velocity[i] - before_body.velocity[i])
                        - external[i];
                    dl[i] += new_l[i] - old_l[i] - applied_l[i];
                }
            }
            for i in 0..3 {
                close(dp[i] + report.planes[1].impulse[i], 0.0, 3e-12);
                close(
                    dl[i] + report.planes[1].angular_impulse[i]
                        - report.coupling.periodic.angular_momentum_change[i],
                    0.0,
                    1e-11,
                );
            }
            close(
                energy(spheres[0]) - energy(old_sphere) + energy(new_carrier)
                    - energy(old_carrier)
                    - work(old_sphere, spheres[0], load)
                    - work(old_carrier, new_carrier, wall_load)
                    - report.coupling.bodies[0].hydrodynamic.body_work
                    + report.planes[1].work
                    + report.coupling.lubrication.dissipated_energy,
                0.0,
                5e-17,
            );
        }
    }
}

#[test]
fn mixed_entry_preserves_empty_prescribed_and_finite_native_results() {
    let collision = CollisionModel3::Bgk { tau: 0.83 };
    let g = grid(collision);
    let fraction = fractions(&g);
    for kind in 0..3 {
        let mut native_grid = grid(collision);
        let mut mixed_grid = grid(collision);
        let mut native_bodies = [sphere()];
        let mut mixed_bodies = native_bodies;
        let mut mixed_surfaces = match kind {
            0 => Vec::new(),
            1 => vec![LubricationSurface3::Prescribed(prescribed())],
            _ => vec![LubricationSurface3::Finite(carrier())],
        };
        let loads = vec![[0.0; 6]; mixed_surfaces.len()];
        let mixed = mixed_grid
            .step_with_lubricated_partial_saturation_surfaces(
                [&fraction],
                [0.7],
                &mut mixed_bodies,
                [[0.0; 6]],
                model(),
                &mut mixed_surfaces,
                &loads,
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        if kind == 0 {
            let native = native_grid
                .step_with_lubricated_partial_saturation_bodies(
                    [&fraction],
                    [0.7],
                    &mut native_bodies,
                    [[0.0; 6]],
                    model(),
                    ORIGIN,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(mixed.coupling, native);
            assert!(mixed.planes.is_empty());
        } else if kind == 1 {
            let native = native_grid
                .step_with_lubricated_partial_saturation_planes(
                    [&fraction],
                    [0.7],
                    &mut native_bodies,
                    [[0.0; 6]],
                    model(),
                    &[prescribed()],
                    ORIGIN,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(mixed, native);
        } else {
            let mut planes = [carrier()];
            let native = native_grid
                .step_with_lubricated_partial_saturation_finite_planes(
                    [&fraction],
                    [0.7],
                    &mut native_bodies,
                    [[0.0; 6]],
                    model(),
                    &mut planes,
                    [[0.0; 6]],
                    ORIGIN,
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(mixed, native);
            assert_eq!(finite(mixed_surfaces[0]), planes[0]);
        }
        assert_eq!(mixed_bodies, native_bodies);
        assert_eq!(snapshot(&mixed_grid), snapshot(&native_grid));
    }
}

#[test]
fn mixed_native_failures_preserve_all_accepted_state_and_retry_matches_fresh() {
    let collision = CollisionModel3::Bgk { tau: 0.83 };
    let mut g = grid(collision);
    let mut fresh = grid(collision);
    let fraction = fractions(&g);
    let before = snapshot(&g);
    let initial_surfaces = surfaces();
    let mut mixed_surfaces = initial_surfaces;
    let mut bodies = [sphere()];
    for (loads, cancelled) in [
        (vec![[0.0; 6]], false),
        (vec![[0.0; 6], [1e-20, 0.0, 0.0, 0.0, 0.0, 0.0]], false),
        (vec![[f64::NAN; 6], [0.0; 6]], false),
        (vec![[1000.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0; 6]], false),
        (vec![[0.0; 6]; 2], true),
    ] {
        let gate = CancelGate::new();
        if cancelled {
            gate.request();
        }
        assert!(
            g.step_with_lubricated_partial_saturation_surfaces(
                [&fraction],
                [0.7],
                &mut bodies,
                [[0.0; 6]],
                model(),
                &mut mixed_surfaces,
                &loads,
                ORIGIN,
                &gate,
            )
            .is_err()
        );
        assert_eq!(bodies, [sphere()]);
        assert_eq!(mixed_surfaces, initial_surfaces);
        assert_eq!(snapshot(&g), before);
    }
    let retry = g
        .step_with_lubricated_partial_saturation_surfaces(
            [&fraction],
            [0.7],
            &mut bodies,
            [[0.0; 6]],
            model(),
            &mut mixed_surfaces,
            &[[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let mut fresh_bodies = [sphere()];
    let mut fresh_surfaces = initial_surfaces;
    let direct = fresh
        .step_with_lubricated_partial_saturation_surfaces(
            [&fraction],
            [0.7],
            &mut fresh_bodies,
            [[0.0; 6]],
            model(),
            &mut fresh_surfaces,
            &[[0.0; 6]; 2],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(retry, direct);
    assert_eq!(bodies, fresh_bodies);
    assert_eq!(mixed_surfaces, fresh_surfaces);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
