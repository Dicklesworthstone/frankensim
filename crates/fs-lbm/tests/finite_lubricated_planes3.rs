//! G1 finite-wall mobility; G3 objective wrench/work; G4 joint rollback.
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, CollisionModel3, E3,
    FiniteLubricationPlane3, LubricatedPlanesExchange3, PartialSaturationBody3, Q3,
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
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn body(center: [f64; 3], mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center,
        mass,
        inertia_world: [
            [0.12, 0.009, -0.004],
            [0.009, 0.17, 0.006],
            [-0.004, 0.006, 0.21],
        ]
        .map(|row| row.map(|x| x * mass)),
        velocity,
        angular_velocity: [0.0; 3],
    }
}
fn grid(model: CollisionModel3) -> BoundaryGrid3 {
    BoundaryGrid3::with_collision_model(8, 8, 8, model, [0.0; 3], BoundarySpec3::periodic())
}
fn kinetic(b: PartialSaturationBody3) -> f64 {
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                b.angular_velocity,
                b.inertia_world.map(|row| dot(row, b.angular_velocity)),
            )
}
fn moment(b: PartialSaturationBody3, origin: [f64; 3]) -> [f64; 3] {
    let orbital = cross(difference(b.center, origin), b.velocity.map(|v| b.mass * v));
    let spin = b.inertia_world.map(|row| dot(row, b.angular_velocity));
    core::array::from_fn(|i| orbital[i] + spin[i])
}
fn external_moment(center: [f64; 3], impulse: [f64; 6], origin: [f64; 3]) -> [f64; 3] {
    let orbital = cross(
        difference(center, origin),
        core::array::from_fn(|i| impulse[i]),
    );
    core::array::from_fn(|i| orbital[i] + impulse[i + 3])
}
fn external_work(
    before: PartialSaturationBody3,
    after: PartialSaturationBody3,
    j: [f64; 6],
) -> f64 {
    (0..3)
        .map(|i| {
            0.5 * j[i] * (before.velocity[i] + after.velocity[i])
                + 0.5 * j[i + 3] * (before.angular_velocity[i] + after.angular_velocity[i])
        })
        .sum()
}
fn frozen<const N: usize, const P: usize>(
    bodies: &mut [PartialSaturationBody3; N],
    external: [[f64; 6]; N],
    radii: [f64; N],
    model: SphereLubrication3,
    planes: &mut [FiniteLubricationPlane3; P],
    plane_external: [[f64; 6]; P],
    origin: [f64; 3],
) -> LubricatedPlanesExchange3<N> {
    // Zero occupancy isolates the common film solve. The native populated
    // body/fluid exchange is independently exercised below.
    let fractions = core::array::from_fn::<_, N, _>(|_| vec![0.0; 512]);
    grid(CollisionModel3::Bgk { tau: 0.83 })
        .step_with_lubricated_partial_saturation_finite_planes(
            fractions.each_ref().map(Vec::as_slice),
            radii,
            bodies,
            external,
            model,
            planes,
            plane_external,
            origin,
            &CancelGate::new(),
        )
        .unwrap()
}

#[test]
fn normal_film_uses_both_inertias_in_one_midpoint_response() {
    let initial = body([3.0, 4.0, 4.0], 2.0, [-0.01, 0.0, 0.0]);
    let wall = body([2.25, 4.0, 4.0], 3.0, [0.015, 0.0, 0.0]);
    let mut spheres = [initial];
    let mut planes =
        [
            FiniteLubricationPlane3::new([2.25, 0.0, 0.0], [1.0, 0.0, 0.0], wall, 0.75, 0.02)
                .unwrap(),
        ];
    let sphere_load = [0.004, 0.0, 0.0, 0.0, 0.0, 0.0];
    let plane_load = [-0.003, 0.0, 0.0, 0.0, 0.0, 0.0];
    let zeta = 6.0 * core::f64::consts::PI * 0.11 * 0.5_f64.powi(2) * (1.0 / 0.25 - 1.0 / 0.75);
    let free_slip = initial.velocity[0] - wall.velocity[0]
        + 0.5 * (sphere_load[0] / initial.mass - plane_load[0] / wall.mass);
    let impulse = -zeta * free_slip / (1.0 + 0.5 * zeta * (1.0 / initial.mass + 1.0 / wall.mass));
    let report = frozen(
        &mut spheres,
        [sphere_load],
        [0.5],
        SphereLubrication3::new(0.11, 0.2, 0.01, 1).unwrap(),
        &mut planes,
        [plane_load],
        ORIGIN,
    );
    let after = planes[0].body();
    close(
        spheres[0].velocity[0],
        initial.velocity[0] + (sphere_load[0] + impulse) / initial.mass,
        2e-17,
    );
    close(
        after.velocity[0],
        wall.velocity[0] + (plane_load[0] - impulse) / wall.mass,
        2e-17,
    );
    close(report.planes[0].impulse[0], -impulse, 3e-17);
    assert!(
        (impulse + zeta * free_slip / (1.0 + 0.5 * zeta / initial.mass)).abs() > 1e-4,
        "a prescribed-wall solve followed by a reaction kick has the wrong impulse"
    );
    let delta = kinetic(spheres[0]) - kinetic(initial) + kinetic(after) - kinetic(wall);
    let work =
        external_work(initial, spheres[0], sphere_load) + external_work(wall, after, plane_load);
    close(
        delta - work + report.coupling.lubrication.dissipated_energy,
        0.0,
        3e-18,
    );
    let response = report.planes[0].finite_body.unwrap();
    assert_eq!(response.body_after, after);
    for r in response.momentum_residual {
        close(r, 0.0, 3e-17);
    }
    close(response.work_residual, 0.0, 3e-18);
    assert_eq!(planes[0].point(), [2.25, 0.0, 0.0]);
    assert_eq!(after.center, wall.center);
    assert_eq!(after.inertia_world, wall.inertia_world);
}

fn shared_plane(wall: PartialSaturationBody3) -> FiniteLubricationPlane3 {
    FiniteLubricationPlane3::new([2.495, 0.0, 0.0], [1.0, 0.0, 0.0], wall, 0.012, 0.001)
        .unwrap()
        .with_tangential(0.01, 0.001)
        .unwrap()
}
fn shared_model() -> SphereLubrication3 {
    SphereLubrication3::new(0.08, 0.4, 0.001, 2)
        .unwrap()
        .with_resistance_mode_limit(10)
        .unwrap()
}

#[test]
fn shared_finite_plane_retains_spin_world_momentum_and_work_at_shifted_origins() {
    let mut initial = [
        body([3.0, 2.9, 4.0], 2.0, [-0.012, 0.007, -0.003]),
        body([3.0, 5.1, 4.2], 3.0, [0.008, -0.009, 0.002]),
    ];
    initial[0].angular_velocity = [0.003, -0.004, 0.002];
    initial[1].angular_velocity = [-0.002, 0.001, -0.005];
    let mut wall = body([2.2, 4.0, 3.7], 4.0, [-0.001, 0.003, 0.002]);
    wall.angular_velocity = [-0.003, 0.004, 0.005];
    let loads = [
        [0.0004, -0.0002, 0.0001, 0.00003, -0.00001, 0.00002],
        [-0.0003, 0.0001, -0.0002, -0.00002, 0.00004, -0.00001],
    ];
    let wall_load = [0.0002, -0.0001, 0.0003, 0.00001, 0.00002, -0.00003];
    let mut accepted = None;
    for origin in [ORIGIN, [8.0, -9.0, 6.0]] {
        let mut spheres = initial;
        let mut planes = [shared_plane(wall)];
        let report = frozen(
            &mut spheres,
            loads,
            [0.5; 2],
            shared_model(),
            &mut planes,
            [wall_load],
            origin,
        );
        let after = planes[0].body();
        let film = report.coupling.lubrication;
        assert_eq!(film.active_plane_pairs, 2);
        assert_eq!(film.active_resistance_modes, 10);
        assert!(film.dissipated_energy > 0.0);
        assert!(
            report.planes[0]
                .finite_body
                .unwrap()
                .spin_impulse
                .iter()
                .any(|v| v.abs() > 1e-5)
        );
        let mut dp = [0.0; 3];
        let mut dl = [0.0; 3];
        let mut energy = 0.0;
        for (before, after, load) in [
            (initial[0], spheres[0], loads[0]),
            (initial[1], spheres[1], loads[1]),
            (wall, after, wall_load),
        ] {
            let old_l = moment(before, origin);
            let new_l = moment(after, origin);
            let imposed_l = external_moment(before.center, load, origin);
            for i in 0..3 {
                dp[i] += before.mass * (after.velocity[i] - before.velocity[i]) - load[i];
                dl[i] += new_l[i] - old_l[i] - imposed_l[i];
            }
            energy += kinetic(after) - kinetic(before) - external_work(before, after, load);
        }
        for v in dp.into_iter().chain(dl) {
            close(v, 0.0, 4e-15);
        }
        close(energy + film.dissipated_energy, 0.0, 8e-18);
        close(film.work_residual, 0.0, 8e-18);
        close(film.constitutive_residual, 0.0, 3e-15);
        for i in 0..3 {
            close(
                film.body_angular_impulses[0][i]
                    + film.body_angular_impulses[1][i]
                    + film.plane_angular_impulse[i],
                0.0,
                3e-15,
            );
        }
        if let Some((old_spheres, old_plane)) = accepted {
            assert_eq!(spheres, old_spheres);
            assert_eq!(after, old_plane);
        }
        accepted = Some((spheres, after));
    }
}

#[test]
fn common_rigid_rotation_is_a_null_mode_of_the_finite_plane_films() {
    let omega = [0.003, -0.002, 0.004];
    let translation = [0.005, -0.003, 0.001];
    let mut initial = [
        body([3.0, 2.9, 4.0], 2.0, [0.0; 3]),
        body([3.0, 5.1, 4.2], 3.0, [0.0; 3]),
    ];
    let mut wall = body([2.2, 4.0, 3.7], 4.0, [0.0; 3]);
    for b in initial.iter_mut().chain(core::iter::once(&mut wall)) {
        b.angular_velocity = omega;
        let spin = cross(omega, difference(b.center, [1.0, 3.0, 2.0]));
        b.velocity = core::array::from_fn(|i| translation[i] + spin[i]);
    }
    let mut spheres = initial;
    let mut planes = [shared_plane(wall)];
    let report = frozen(
        &mut spheres,
        [[0.0; 6]; 2],
        [0.5; 2],
        shared_model(),
        &mut planes,
        [[0.0; 6]],
        ORIGIN,
    );
    for (before, after) in initial
        .into_iter()
        .zip(spheres)
        .chain(core::iter::once((wall, planes[0].body())))
    {
        for i in 0..3 {
            close(before.velocity[i], after.velocity[i], 1e-16);
            close(before.angular_velocity[i], after.angular_velocity[i], 1e-16);
        }
    }
    close(report.coupling.lubrication.dissipated_energy, 0.0, 1e-30);
}

fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..512)
        .map(|i| g.populations(i % 8, (i / 8) % 8, i / 64))
        .collect()
}

#[test]
fn finite_recoil_and_forced_native_fluid_close_one_shared_tick() {
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
        let mut initial = body([3.5, 4.5, 4.5], 3.0, [-0.008, 0.004, -0.002]);
        initial.angular_velocity = [0.003, -0.004, 0.002];
        let wall = body([2.6, 4.1, 4.8], 5.0, [0.002, -0.001, 0.003]);
        let mut spheres = [initial];
        let mut planes =
            [
                FiniteLubricationPlane3::new([2.79, 0.0, 0.0], [1.0, 0.0, 0.0], wall, 0.03, 0.001)
                    .unwrap()
                    .with_tangential(0.013, 0.001)
                    .unwrap(),
            ];
        let mut sampler = BodyFrameSampler3::new(&g, 4).unwrap();
        sampler
            .sample(
                initial.center,
                So3::identity(),
                &BodyFrameEllipsoid3::sphere(0.7).unwrap(),
                &CancelGate::new(),
            )
            .unwrap();
        let before = snapshot(&g);
        let report = g
            .step_with_lubricated_partial_saturation_finite_planes(
                [sampler.fractions()],
                [0.7],
                &mut spheres,
                [[0.0; 6]],
                SphereLubrication3::new(0.11, 0.1, 0.001, 1).unwrap(),
                &mut planes,
                [[0.0; 6]],
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        assert!(report.coupling.bodies[0].hydrodynamic.coupled_cells > 0);
        let after = snapshot(&g);
        let mut dp = [0.0; 3];
        let mut dl = [0.0; 3];
        for i in 0..512 {
            let cell = [i % 8, (i / 8) % 8, i / 64];
            let mut local = [0.0; 3];
            for q in 0..Q3 {
                let delta = after[i][q] - before[i][q];
                for (a, c) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                    local[a] += f64::from(c) * delta;
                }
            }
            if cell.into_iter().all(|j| (2..6).contains(&j)) {
                for a in 0..3 {
                    local[a] -= force[a];
                }
            }
            let torque = cross(
                core::array::from_fn(|a| cell[a] as f64 + 0.5 - ORIGIN[a]),
                local,
            );
            for a in 0..3 {
                dp[a] += local[a];
                dl[a] += torque[a];
            }
        }
        for (before, after) in [(initial, spheres[0]), (wall, planes[0].body())] {
            let old_l = moment(before, ORIGIN);
            let new_l = moment(after, ORIGIN);
            for a in 0..3 {
                dp[a] += before.mass * (after.velocity[a] - before.velocity[a]);
                dl[a] += new_l[a] - old_l[a];
            }
        }
        for a in 0..3 {
            close(dp[a], 0.0, 3e-12);
            close(
                dl[a] - report.coupling.periodic.angular_momentum_change[a],
                0.0,
                1e-11,
            );
        }
        let kinetic_delta =
            kinetic(spheres[0]) - kinetic(initial) + kinetic(planes[0].body()) - kinetic(wall);
        close(
            kinetic_delta + report.coupling.lubrication.dissipated_energy
                - report.coupling.bodies[0].hydrodynamic.body_work,
            0.0,
            4e-17,
        );
    }
}

#[test]
fn inactive_finite_surface_accepts_external_impulse_and_refusal_is_atomic() {
    let initial = body([3.0, 4.0, 4.0], 2.0, [0.0; 3]);
    let wall = body([0.0, 4.0, 4.0], 3.0, [0.0; 3]);
    let plane = FiniteLubricationPlane3::new([0.0; 3], [1.0, 0.0, 0.0], wall, 0.2, 0.01).unwrap();
    let load = [0.003, -0.002, 0.001, 0.0004, -0.0002, 0.0003];
    let mut spheres = [initial];
    let mut planes = [plane];
    let model = SphereLubrication3::new(0.1, 0.2, 0.01, 1).unwrap();
    let report = frozen(
        &mut spheres,
        [[0.0; 6]],
        [0.5],
        model,
        &mut planes,
        [load],
        ORIGIN,
    );
    assert_eq!(report.coupling.lubrication.active_resistance_modes, 0);
    let after = planes[0].body();
    let dl = wall.inertia_world.map(|row| {
        dot(
            row,
            difference(after.angular_velocity, wall.angular_velocity),
        )
    });
    for i in 0..3 {
        close(wall.mass * after.velocity[i], load[i], 2e-18);
        close(dl[i], load[i + 3], 2e-18);
    }
    close(
        kinetic(after) - kinetic(wall) - external_work(wall, after, load),
        0.0,
        2e-20,
    );
    let mut g = grid(CollisionModel3::Bgk { tau: 0.83 });
    let before = snapshot(&g);
    let fractions = vec![0.0; 512];
    for cancelled in [false, true] {
        let mut spheres = [initial];
        let mut planes = [plane];
        let gate = CancelGate::new();
        if cancelled {
            gate.request();
        }
        let failed = g.step_with_lubricated_partial_saturation_finite_planes(
            [&fractions],
            [0.5],
            &mut spheres,
            [[0.0; 6]],
            model,
            &mut planes,
            [[1000.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
            ORIGIN,
            &gate,
        );
        assert!(failed.is_err());
        assert_eq!(spheres, [initial]);
        assert_eq!(planes, [plane]);
        assert_eq!(snapshot(&g), before);
    }
}
