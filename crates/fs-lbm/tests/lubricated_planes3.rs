//! G1 normal wall asymptotics, G3 coupled reactions and G4 atomic geometry.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, E3, FaceBoundary3, LubricatedPlanesExchange3,
    LubricationPlane3, MovingPartialSaturationBody3, PartialSaturationBody3,
    PartialSaturationError3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}: {:.3e}",
        a - b
    );
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
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn grid(walls: bool) -> BoundaryGrid3 {
    BoundaryGrid3::new(
        8,
        8,
        8,
        0.83,
        [0.0; 3],
        if walls {
            BoundarySpec3::new([FaceBoundary3::Wall { velocity: [0.0; 3] }; 6])
        } else {
            BoundarySpec3::periodic()
        },
    )
}
fn body(x: f64, mass: f64, velocity: [f64; 3]) -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [x, 4.0, 4.0],
        mass,
        inertia_world: [
            [0.006, 0.0004, -0.0002],
            [0.0004, 0.008, 0.0003],
            [-0.0002, 0.0003, 0.009],
        ],
        velocity,
        angular_velocity: [0.002, -0.001, 0.003],
    }
}
fn plane(x: f64, normal: f64, speed: f64, cutoff: f64, roughness: f64) -> LubricationPlane3 {
    LubricationPlane3::new(
        [x, 0.0, 0.0],
        [normal, 0.0, 0.0],
        [speed, 0.0, 0.0],
        cutoff,
        roughness,
    )
    .unwrap()
}
fn zeta(mu: f64, radius: f64, gap: f64, cutoff: f64, roughness: f64) -> f64 {
    6.0 * core::f64::consts::PI
        * mu
        * radius
        * radius
        * (1.0 / gap.max(roughness) - 1.0 / cutoff).max(0.0)
}
fn frozen<const N: usize>(
    mut bodies: [PartialSaturationBody3; N],
    external: [[f64; 6]; N],
    radii: [f64; N],
    model: SphereLubrication3,
    planes: &[LubricationPlane3],
) -> ([PartialSaturationBody3; N], LubricatedPlanesExchange3<N>) {
    // Zero-occupancy isolates the subgrid resistance from native resolved drag.
    // A fully sampled sphere beside actual voxel walls is exercised below.
    let fields = core::array::from_fn::<_, N, _>(|_| vec![0.0; 512]);
    let report = grid(false)
        .step_with_lubricated_partial_saturation_planes(
            fields.each_ref().map(Vec::as_slice),
            radii,
            &mut bodies,
            external,
            model,
            planes,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    (bodies, report)
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                result.push(if g.is_solid(x, y, z) {
                    [0.0; Q3]
                } else {
                    g.populations(x, y, z)
                });
            }
        }
    }
    result
}

#[test]
fn driven_plane_matches_stokes_solution_world_moment_and_source_work() {
    let r = 0.5;
    let resistance = zeta(0.11, r, 0.25, 0.75, 0.02);
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.31, -0.27, 0.19))).unwrap();
    for rotated in [false, true] {
        for boost in [0.0, 0.011] {
            let transform = |x| {
                if rotated {
                    array(rotation.rotate(vector(x)).unwrap())
                } else {
                    x
                }
            };
            let mut initial = body(3.0, 2.0, [0.003 + boost, -0.002, 0.001]);
            initial.center = transform(initial.center);
            initial.velocity = transform(initial.velocity);
            initial.angular_velocity = transform(initial.angular_velocity);
            if rotated {
                let tensor = Mat3 {
                    m: core::array::from_fn(|i| initial.inertia_world[i / 3][i % 3]),
                };
                let q = rotation.matrix();
                let tensor = q.compose(tensor).compose(q.transpose());
                initial.inertia_world = core::array::from_fn(|i| {
                    core::array::from_fn(|j| tensor.m[3 * i.max(j) + i.min(j)])
                });
            }
            let force = transform([0.004, -0.002, 0.001]);
            let n = transform([1.0, 0.0, 0.0]);
            let wall = LubricationPlane3::new(
                transform([2.25, 0.0, 0.0]),
                n,
                transform([0.02 + boost, 0.004, -0.003]),
                0.75,
                0.02,
            )
            .unwrap();
            let lambda =
                -resistance * (0.003 - 0.02 + 0.5 * 0.004 / 2.0) / (1.0 + 0.5 * resistance / 2.0);
            let external = [force[0], force[1], force[2], 0.0, 0.0, 0.0];
            let (after, report) = frozen(
                [initial],
                [external],
                [r],
                SphereLubrication3::new(0.11, 0.2, 0.01, 1).unwrap(),
                &[wall],
            );
            let film = report.coupling.lubrication;
            assert_eq!(film.active_pairs, 0);
            assert_eq!(film.active_plane_pairs, 1);
            assert_eq!(report.planes[0].active_spheres, 1);
            for i in 0..3 {
                close(
                    after[0].velocity[i],
                    initial.velocity[i] + (force[i] + lambda * n[i]) / 2.0,
                    7e-17,
                );
                close(
                    after[0].angular_velocity[i],
                    initial.angular_velocity[i],
                    2e-17,
                );
                close(film.body_impulses[0][i], lambda * n[i], 1e-16);
                close(report.planes[0].impulse[i], -lambda * n[i], 1e-16);
                close(film.body_impulses[0][i] + film.plane_impulse[i], 0.0, 0.0);
                close(
                    film.body_angular_impulses[0][i] + film.plane_angular_impulse[i],
                    0.0,
                    1e-16,
                );
                close(report.coupling.bodies[0].momentum_residual[i], 0.0, 2e-17);
            }
            let midpoint =
                core::array::from_fn(|i| 0.5 * (initial.velocity[i] + after[0].velocity[i]));
            let actual_work = dot(film.body_impulses[0], midpoint);
            let plane_work = dot(report.planes[0].impulse, wall.velocity());
            assert!(actual_work > 0.0, "the moving support supplies real power");
            assert!(plane_work < 0.0);
            close(report.planes[0].work, plane_work, 2e-18);
            close(film.body_work[0], actual_work, 2e-18);
            close(
                actual_work + plane_work + film.dissipated_energy,
                0.0,
                3e-18,
            );
            close(film.work_residual, 0.0, 3e-18);
            close(film.constitutive_residual, 0.0, 7e-17);
            let kinetic =
                dot(after[0].velocity, after[0].velocity) - dot(initial.velocity, initial.velocity);
            close(kinetic, dot(force, midpoint) + actual_work, 3e-18);
        }
    }
}

#[test]
fn sphere_pair_and_two_driven_planes_are_one_simultaneous_system() {
    let initial = [
        body(3.0, 2.0, [-0.012, 0.003, 0.0]),
        body(4.2, 3.0, [0.018, 0.003, 0.0]),
    ];
    let external = [
        [0.005, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.003, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let planes = [
        plane(2.3, 1.0, 0.008, 0.45, 0.03),
        plane(4.9, -1.0, -0.006, 0.7, 0.02),
    ];
    let pair = zeta(0.11, 0.25, 0.2, 0.6, 0.01);
    let left = zeta(0.11, 0.5, 0.2, 0.45, 0.03);
    let right = zeta(0.11, 0.5, 0.2, 0.7, 0.02);
    // Independently assemble the TWO body-velocity equations; production solves
    // THREE film impulses using factored six-velocity native body responses.
    let a = 2.0 + 0.5 * (left + pair);
    let d = 3.0 + 0.5 * (right + pair);
    let b = -0.5 * pair;
    let x = (2.0 - 0.5 * (left + pair)) * initial[0].velocity[0]
        + 0.5 * pair * initial[1].velocity[0]
        + external[0][0]
        + left * 0.008;
    let y = (3.0 - 0.5 * (right + pair)) * initial[1].velocity[0]
        + 0.5 * pair * initial[0].velocity[0]
        + external[1][0]
        - right * 0.006;
    let expected = [
        (d * x - b * y) / (a * d - b * b),
        (a * y - b * x) / (a * d - b * b),
    ];
    for reverse in [false, true] {
        let order = if reverse { [1, 0] } else { [0, 1] };
        let (after, result) = frozen(
            order.map(|i| initial[i]),
            order.map(|i| external[i]),
            [0.5; 2],
            SphereLubrication3::new(0.11, 0.6, 0.01, 3).unwrap(),
            &order.map(|i| planes[i]),
        );
        assert_eq!(result.coupling.lubrication.active_pairs, 1);
        assert_eq!(result.coupling.lubrication.active_plane_pairs, 2);
        for i in 0..2 {
            close(after[i].velocity[0], expected[order[i]], 8e-17);
        }
        let film = result.coupling.lubrication;
        for axis in 0..3 {
            close(
                film.body_impulses[0][axis]
                    + film.body_impulses[1][axis]
                    + film.plane_impulse[axis],
                0.0,
                3e-17,
            );
            close(
                film.body_angular_impulses[0][axis]
                    + film.body_angular_impulses[1][axis]
                    + film.plane_angular_impulse[axis],
                0.0,
                1e-16,
            );
        }
        close(film.work_residual, 0.0, 3e-18);
        assert!(film.dissipated_energy > 0.0);
    }
}

#[test]
fn roughness_cutoff_and_empty_plane_paths_have_the_declared_limits() {
    let model = SphereLubrication3::new(0.11, 0.7, 0.01, 64).unwrap();
    let initial = body(0.75, 1.0, [-0.01, 0.003, 0.0]);
    let r = 0.25;
    let mut previous = None;
    for gap in [0.015625, 0.0625] {
        let wall = plane(initial.center[0] - r - gap, 1.0, 0.0, 0.5, 0.125);
        let (after, report) = frozen([initial], [[0.0; 6]], [r], model, &[wall]);
        if let Some(expected) = previous {
            assert_eq!(after, expected);
        }
        previous = Some(after);
        let resistance = zeta(0.11, r, 0.125, 0.5, 0.125);
        close(
            after[0].velocity[0],
            initial.velocity[0] * (1.0 - 0.5 * resistance) / (1.0 + 0.5 * resistance),
            2e-17,
        );
        assert_eq!(report.coupling.lubrication.active_plane_pairs, 1);
    }
    let (_, at_cutoff) = frozen(
        [initial],
        [[0.0; 6]],
        [r],
        model,
        &[plane(0.0, 1.0, 0.0, 0.5, 0.125)],
    );
    assert_eq!(at_cutoff.coupling.lubrication.active_plane_pairs, 0);
    assert_eq!(at_cutoff.planes[0].impulse, [0.0; 3]);
    let fields = [vec![0.0; 512]; 1];
    let mut a = [initial];
    let mut b = a;
    let mut ga = grid(false);
    let mut gb = grid(false);
    let old = ga
        .step_with_lubricated_partial_saturation_bodies(
            fields.each_ref().map(Vec::as_slice),
            [r],
            &mut a,
            [[0.0; 6]],
            model,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let new = gb
        .step_with_lubricated_partial_saturation_planes(
            fields.each_ref().map(Vec::as_slice),
            [r],
            &mut b,
            [[0.0; 6]],
            model,
            &[],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(old, new.coupling);
    assert!(new.planes.is_empty());
    assert_eq!(snapshot(&ga), snapshot(&gb));
    let mut a = [moving(0), moving(1)];
    let mut b = a;
    let mut ga = grid(false);
    let mut gb = grid(false);
    let mut sa = core::array::from_fn(|_| BodyFrameSampler3::new(&ga, 3).unwrap());
    let mut sb = core::array::from_fn(|_| BodyFrameSampler3::new(&gb, 3).unwrap());
    let old = ga
        .step_with_lubricated_spheres(
            &mut sa,
            [0.5; 2],
            &mut a,
            [[0.0; 6]; 2],
            model,
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let new = gb
        .step_with_lubricated_spheres_and_planes(
            &mut sb,
            [0.5; 2],
            &mut b,
            [[0.0; 6]; 2],
            model,
            &[],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(old, new.motion);
    assert!(new.planes.is_empty());
    assert_eq!(snapshot(&ga), snapshot(&gb));
}

fn moving(a: usize) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [2.7 + 1.15 * a as f64, 4.0, 4.0],
        orientation: So3::identity(),
        mass: 1.3 + 0.2 * a as f64,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                0.06, 0.008, -0.004, 0.008, 0.08, 0.003, -0.004, 0.003, 0.095,
            ],
        })
        .unwrap(),
        velocity: [-0.003 - 0.001 * a as f64, 0.001, -0.001],
        angular_momentum_body: Vec3::new(0.00008, -0.00006, 0.00004),
    }
}
fn momentum(b: MovingPartialSaturationBody3) -> ([f64; 3], [f64; 3]) {
    let p = b.velocity.map(|v| b.mass * v);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    let orbital = cross(core::array::from_fn(|i| b.center[i] - ORIGIN[i]), p);
    (p, core::array::from_fn(|i| spin[i] + orbital[i]))
}
fn fluid_delta(before: &[[f64; Q3]], g: &BoundaryGrid3) -> ([f64; 3], [f64; 3]) {
    let (mut p, mut h) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if g.is_solid(x, y, z) {
                    continue;
                }
                let f = g.populations(x, y, z);
                let mut local = [0.0; 3];
                for q in 0..Q3 {
                    let d = f[q] - before[(z * 8 + y) * 8 + x][q];
                    for (i, e) in [E3[q].0, E3[q].1, E3[q].2].into_iter().enumerate() {
                        local[i] += f64::from(e) * d;
                    }
                }
                let torque = cross(
                    [
                        x as f64 + 0.5 - ORIGIN[0],
                        y as f64 + 0.5 - ORIGIN[1],
                        z as f64 + 0.5 - ORIGIN[2],
                    ],
                    local,
                );
                for i in 0..3 {
                    p[i] += local[i];
                    h[i] += torque[i];
                }
            }
        }
    }
    (p, h)
}

#[test]
fn sampled_full_tensor_spheres_beside_voxel_wall_close_forced_world_balances() {
    let mut g = grid(true);
    // Cell centers x=0.5 and 1.5 are solid; the actual halfway wall is x=2.
    g.voxelize_sdf(|p| p[0] - 2.0);
    let force = [-2e-7, 1e-7, 0.5e-7];
    g.set_force_field(|_| force);
    let (mut forcing, mut force_moment) = ([0.0; 3], [0.0; 3]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if g.is_solid(x, y, z) {
                    continue;
                }
                let torque = cross(
                    [
                        x as f64 + 0.5 - ORIGIN[0],
                        y as f64 + 0.5 - ORIGIN[1],
                        z as f64 + 0.5 - ORIGIN[2],
                    ],
                    force,
                );
                for i in 0..3 {
                    forcing[i] += force[i];
                    force_moment[i] += torque[i];
                }
            }
        }
    }
    let mut bodies = [moving(0), moving(1)];
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&g, 4).unwrap());
    let model = SphereLubrication3::new(0.11, 0.5, 0.02, 3).unwrap();
    let wall = plane(2.0, 1.0, 0.0, 0.6, 0.03);
    let external = [[-0.00001, 0.000004, -0.000002, 0.000001, -0.000002, 0.000001]; 2];
    let mut maximum_defect = 0.0_f64;
    let mut loss = 0.0;
    for _ in 0..24 {
        let old = bodies;
        let before = snapshot(&g);
        let result = g
            .step_with_lubricated_spheres_and_planes(
                &mut samplers,
                [0.5; 2],
                &mut bodies,
                external,
                model,
                &[wall],
                ORIGIN,
                &CancelGate::new(),
            )
            .unwrap();
        let report = result.motion;
        let film = report.coupling.lubrication;
        assert_eq!(film.active_pairs, 1);
        assert_eq!(film.active_plane_pairs, 1);
        assert_eq!(result.planes[0].active_spheres, 1);
        let (mut dp, mut dh) = fluid_delta(&before, &g);
        for a in 0..2 {
            let previous = momentum(old[a]);
            let now = momentum(bodies[a]);
            let torque = cross(
                core::array::from_fn(|i| report.bodies[a].midpoint_center[i] - ORIGIN[i]),
                external[a][..3].try_into().unwrap(),
            );
            for i in 0..3 {
                dp[i] += now.0[i] - previous.0[i] - external[a][i];
                dh[i] += now.1[i] - previous.1[i] - external[a][i + 3] - torque[i];
                close(report.bodies[a].momentum_residual[i], 0.0, 2e-15);
                close(report.bodies[a].momentum_residual[i + 3], 0.0, 2e-15);
            }
            close(report.bodies[a].work_residual, 0.0, 3e-17);
            assert!(bodies[a].center[0] - 0.5 > 2.0);
        }
        let h = report.bodies[0].midpoint_center[0] - 2.0 - 0.5;
        let vmid = 0.5 * (old[0].velocity[0] + bodies[0].velocity[0]);
        let resistance = zeta(0.11, 0.5, h, 0.6, 0.03);
        close(result.planes[0].impulse[0], resistance * vmid, 2e-17);
        close(film.plane_work, 0.0, 0.0);
        close(film.work_residual, 0.0, 3e-18);
        loss += film.dissipated_energy;
        for i in 0..3 {
            let p =
                dp[i] + report.coupling.walls.wall_impulse[i] + film.plane_impulse[i] - forcing[i];
            let h = dh[i]
                + report.coupling.walls.wall_angular_impulse[i]
                + film.plane_angular_impulse[i]
                - force_moment[i];
            close(p, 0.0, 3e-12);
            close(h, 0.0, 9e-12);
            maximum_defect = maximum_defect.max(p.abs()).max(h.abs());
        }
    }
    assert!(loss > 0.0);
    eprintln!(
        "24 sampled wall-film ticks: maximum world momentum defect={maximum_defect:e}, modeled loss={loss:e}"
    );
}

#[test]
fn moving_plane_geometry_uses_both_half_drifts_and_late_failure_is_atomic() {
    let mut g = grid(false);
    let initial_fluid = snapshot(&g);
    let mut initial = moving(0);
    initial.center = [3.0, 4.0, 4.0];
    initial.velocity = [0.0; 3];
    initial.angular_momentum_body = Vec3::new(0.0, 0.0, 0.0);
    let mut bodies = [initial];
    let mut samplers = [BodyFrameSampler3::new(&g, 3).unwrap()];
    let model = SphereLubrication3::new(1e-10, 0.5, 0.01, 1).unwrap();
    // Start gap .02 and midpoint gap .005 are legal, but the prescribed wall
    // crosses the almost-stationary sphere during the second half-drift.
    let wall = plane(2.48, 1.0, 0.03, 0.5, 0.01);
    let refused = g.step_with_lubricated_spheres_and_planes(
        &mut samplers,
        [0.5],
        &mut bodies,
        [[0.0; 6]],
        model,
        &[wall],
        ORIGIN,
        &CancelGate::new(),
    );
    assert!(matches!(
        refused,
        Err(PartialSaturationError3::Invalid {
            reason: "moving sphere meets prescribed lubrication plane"
        })
    ));
    assert_eq!(bodies, [initial]);
    assert_eq!(snapshot(&g), initial_fluid);
    let wall = plane(2.48, 1.0, 0.001, 0.5, 0.01);
    let accepted = g
        .step_with_lubricated_spheres_and_planes(
            &mut samplers,
            [0.5],
            &mut bodies,
            [[0.0; 6]],
            model,
            &[wall],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    let mut fresh = grid(false);
    let mut expected = [initial];
    let mut fresh_sampler = [BodyFrameSampler3::new(&fresh, 3).unwrap()];
    let repeated = fresh
        .step_with_lubricated_spheres_and_planes(
            &mut fresh_sampler,
            [0.5],
            &mut expected,
            [[0.0; 6]],
            model,
            &[wall],
            ORIGIN,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(bodies, expected);
    assert_eq!(snapshot(&g), snapshot(&fresh));
    // Reaction uses the translated midpoint, not the initial .02 gap.
    let resistance = zeta(1e-10, 0.5, 0.0195, 0.5, 0.01);
    let slip = 0.5 * (initial.velocity[0] + bodies[0].velocity[0]) - 0.001;
    close(accepted.planes[0].impulse[0], resistance * slip, 2e-22);
    close(wall.translated(1.0).unwrap().point()[0], 2.481, 1e-15);
}

#[test]
fn plane_admission_combined_pair_caps_and_cancellation_preserve_state() {
    for bad in [f64::NAN, f64::INFINITY, -0.1, 0.0] {
        assert!(LubricationPlane3::new([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], bad, 0.01).is_err());
    }
    assert!(LubricationPlane3::new([0.0; 3], [0.0; 3], [0.0; 3], 0.5, 0.01).is_err());
    assert!(LubricationPlane3::new([0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 0.01).is_err());
    let mut g = grid(false);
    let old_fluid = snapshot(&g);
    let initial = [
        body(3.0, 2.0, [-0.01, 0.0, 0.0]),
        body(4.2, 3.0, [0.01, 0.0, 0.0]),
    ];
    let mut bodies = initial;
    let fields = [vec![0.0; 512], vec![0.0; 512]];
    let model = SphereLubrication3::new(0.11, 0.6, 0.01, 2).unwrap();
    let left = plane(2.3, 1.0, 0.0, 0.5, 0.01);
    let right = plane(4.9, -1.0, 0.0, 0.5, 0.01);
    for planes in [
        vec![left, right],
        vec![left, left],
        vec![plane(2.51, 1.0, 0.0, 0.001, 0.001)],
    ] {
        assert!(
            g.step_with_lubricated_partial_saturation_planes(
                fields.each_ref().map(Vec::as_slice),
                [0.5; 2],
                &mut bodies,
                [[0.0; 6]; 2],
                model,
                &planes,
                ORIGIN,
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(bodies, initial);
        assert_eq!(snapshot(&g), old_fluid);
    }
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        g.step_with_lubricated_partial_saturation_planes(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            model,
            &[left],
            ORIGIN,
            &gate
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&g), old_fluid);
    let planes = (0..64)
        .map(|i| plane(-1.0 - i as f64, 1.0, 0.0, 0.5, 0.01))
        .collect::<Vec<_>>();
    let (_, report) = frozen([initial[0]], [[0.0; 6]], [0.5], model, &planes);
    assert_eq!(report.planes.len(), 64);
    assert_eq!(report.coupling.lubrication.active_plane_pairs, 0);
    let mut too_many = planes;
    too_many.push(plane(-70.0, 1.0, 0.0, 0.5, 0.01));
    assert!(
        g.step_with_lubricated_partial_saturation_planes(
            fields.each_ref().map(Vec::as_slice),
            [0.5; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            model,
            &too_many,
            ORIGIN,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&g), old_fluid);
}
