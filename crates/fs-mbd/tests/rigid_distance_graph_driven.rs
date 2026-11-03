use fs_mbd::rigid_rattle::distance::graph::{
    DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenDistanceGraphStep, advance_distance_graph, advance_driven_distance_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    near(a.x, b.x, tol);
    near(a.y, b.y, tol);
    near(a.z, b.z, tol);
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn controls() -> DistanceGraphControls {
    DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 10.,
            maximum_rotation_increment_rad: 0.25,
            maximum_position_iterations: 40,
            maximum_backtracks: 16,
        },
        position_tolerance_m: 2e-13,
        velocity_tolerance_mps: 2e-13,
        minimum_scaled_pivot: 1e-10,
        maximum_bodies: 64,
        maximum_constraints: 32,
        maximum_trial_evaluations: 160,
        maximum_live_bytes: 2_000_000,
    }
}
fn spherical(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(0.4, 0.4, 0.4).scale(mass)).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn center_link(a: usize, b: usize, length: f64) -> DistanceLink {
    DistanceLink {
        bodies: [a, b],
        points_body_m: [Vec3::ZERO; 2],
        length_m: length,
    }
}
fn run(
    s: &[RigidBodyState],
    p: &[MassProperties],
    links: &[DistanceLink],
    drives: &[DistanceDrive],
    h: f64,
) -> DrivenDistanceGraphStep {
    advance_driven_distance_graph(
        s,
        p,
        links,
        drives,
        free(),
        &vec![Wrench::ZERO; s.len()],
        h,
        controls(),
        || false,
    )
    .unwrap()
}
fn totals(s: &[RigidBodyState], p: &[MassProperties]) -> (Vec3, Vec3, f64) {
    s.iter().zip(p).fold(
        (Vec3::ZERO, Vec3::ZERO, 0.),
        |(linear, angular, energy), (&s, &p)| {
            let d = free().diagnostics(s, p).unwrap();
            (
                linear.add(d.linear_momentum_world),
                angular.add(d.angular_momentum_world),
                energy + d.mechanical_energy,
            )
        },
    )
}
fn state_near(a: RigidBodyState, b: RigidBodyState, tol: f64) {
    vector(a.pose().position_world(), b.pose().position_world(), tol);
    vector(a.linear_momentum_world(), b.linear_momentum_world(), tol);
    vector(a.angular_momentum_body(), b.angular_momentum_body(), tol);
    for e in [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ] {
        vector(
            a.pose().orientation().rotate_body_to_world(e),
            b.pose().orientation().rotate_body_to_world(e),
            tol,
        );
    }
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.),
        s.angular_momentum_body().scale(-1.),
    )
    .unwrap()
}
fn chain_rates(r: [f64; 2]) -> [f64; 3] {
    let first = -(2. * r[0] + 3. * (r[0] + r[1])) / 6.;
    [first, first + r[0], first + r[0] + r[1]]
}

#[test]
fn g1_driven_chain_shares_recoil_and_signed_actuator_work_analytically() {
    let p = [spherical(1.), spherical(2.), spherical(3.)];
    let n = Vec3::new(0.6, 0.8, 0.);
    let boost = Vec3::new(0.2, -0.3, 0.4);
    let links = [center_link(0, 1, 1.), center_link(1, 2, 2.)];
    let speeds = [0.3, -0.2];
    let velocities = chain_rates(speeds);
    let s = std::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            n.scale([0., 1., 3.][i]),
            UnitQuaternion::IDENTITY,
            boost.add(n.scale(velocities[i])),
            Vec3::ZERO,
        )
    });
    let before = totals(&s, &p);
    let h = 0.2;
    for sign in [1., -1.] {
        let accelerations = [0.4 * sign, 0.1 * sign];
        let body_accelerations = chain_rates(accelerations);
        let drives = std::array::from_fn::<_, 2, _>(|i| DistanceDrive {
            range_rate_mps: speeds[i],
            range_acceleration_mps2: accelerations[i],
        });
        let out = run(&s, &p, &links, &drives, h);
        for i in 0..3 {
            vector(
                out.step.states_after[i].pose().position_world(),
                s[i].pose()
                    .position_world()
                    .add(boost.add(n.scale(velocities[i])).scale(h))
                    .add(n.scale(0.5 * h * h * body_accelerations[i])),
                3e-14,
            );
            vector(
                out.step.states_after[i]
                    .center_of_mass_velocity_world(p[i])
                    .unwrap(),
                boost.add(n.scale(velocities[i] + h * body_accelerations[i])),
                3e-14,
            );
        }
        // Independent finite-mass force balance: each endpoint receives a
        // half-step of its constant acceleration, with both middle-body terms.
        let impulses = [
            -0.5 * h * body_accelerations[0],
            0.5 * h * 3. * body_accelerations[2],
        ];
        let mut exact_work = 0.;
        for i in 0..2 {
            for events in [&out.step.position_impulses, &out.step.velocity_impulses] {
                vector(events[i].second.impulse_world, n.scale(impulses[i]), 3e-14);
            }
            let work = impulses[i] * (2. * speeds[i] + h * accelerations[i]);
            exact_work += work;
            near(out.links[i].actuator_work_j, work, 3e-14);
            near(
                out.links[i].link_after.length_m,
                links[i].length_m + speeds[i] * h + 0.5 * accelerations[i] * h * h,
                2e-15,
            );
            near(
                out.links[i].drive_after.range_rate_mps,
                speeds[i] + accelerations[i] * h,
                2e-16,
            );
        }
        let after = totals(&out.step.states_after, &p);
        vector(after.0, before.0, 2e-15);
        vector(after.1, before.1, 3e-14);
        near(after.2 - before.2, exact_work, 3e-14);
        near(out.actuator_work_j, exact_work, 3e-14);
        near(out.constraint_work_defect_j, 0., 3e-14);
        assert_eq!(out.actuator_work_j > 0., sign > 0.);
    }
}

#[test]
fn g0_zero_commands_and_unlinked_bodies_keep_the_exact_native_path() {
    let p = [spherical(1.), spherical(2.), spherical(3.)];
    let s = std::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            Vec3::new(i as f64, 0., 0.),
            UnitQuaternion::IDENTITY,
            Vec3::new(0., 0.1, 0.2),
            Vec3::ZERO,
        )
    });
    let links = [center_link(0, 1, 1.), center_link(1, 2, 1.)];
    let forces = [7., -2., 1.].map(|f| Wrench {
        force_world: Vec3::new(f, 0., 0.),
        torque_body: Vec3::ZERO,
    });
    for count in 0..=2 {
        let old = advance_distance_graph(
            &s,
            &p,
            &links[..count],
            free(),
            &forces,
            0.1,
            controls(),
            || false,
        )
        .unwrap();
        let driven = advance_driven_distance_graph(
            &s,
            &p,
            &links[..count],
            &vec![DistanceDrive::default(); count],
            free(),
            &forces,
            0.1,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(old, driven.step);
        near(driven.actuator_work_j, 0., 0.);
        near(
            driven.constraint_work_defect_j,
            old.constraint_work_j,
            2e-16,
        );
    }
}

#[test]
fn g0_driven_spin_survives_common_boosts_before_rate_and_work_pairing() {
    let properties = [spherical(1.), spherical(2.)];
    let link = DistanceLink {
        bodies: [0, 1],
        points_body_m: [Vec3::ZERO, Vec3::new(0., 1., 0.)],
        length_m: 2_f64.sqrt(),
    };
    let scene = |boost| {
        [
            state(
                properties[0],
                Vec3::ZERO,
                UnitQuaternion::IDENTITY,
                Vec3::new(boost, 0., 0.),
                Vec3::ZERO,
            ),
            state(
                properties[1],
                Vec3::new(1., 0., 0.),
                UnitQuaternion::IDENTITY,
                Vec3::new(boost, 0., 0.),
                Vec3::new(0., 0., 1.),
            ),
        ]
    };
    let hidden = scene(2_f64.powi(54));
    let first = hidden[0]
        .point_kinematics(properties[0], link.points_body_m[0])
        .unwrap();
    let second = hidden[1]
        .point_kinematics(properties[1], link.points_body_m[1])
        .unwrap();
    assert_eq!(first.point_velocity_world, second.point_velocity_world);
    near(
        second.angular_velocity_world.cross(second.arm_world).x,
        -1.,
        0.,
    );
    assert!(matches!(
        advance_driven_distance_graph(
            &hidden,
            &properties,
            &[link],
            &[DistanceDrive::default()],
            free(),
            &[Wrench::ZERO; 2],
            0.001,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { constraint: 0 })
    ));

    // At boosts where the actual reaction increments remain representable,
    // compare motion and source/defect work with the unboosted native result.
    let command = DistanceDrive {
        range_rate_mps: -1. / 2_f64.sqrt(),
        range_acceleration_mps2: 0.15,
    };
    let h = 0.001;
    let baseline = run(&scene(0.), &properties, &[link], &[command], h);
    for boost in [8., 64.] {
        let out = run(&scene(boost), &properties, &[link], &[command], h);
        for i in 0..2 {
            let expected = transform(
                baseline.step.states_after[i],
                properties[i],
                UnitQuaternion::IDENTITY,
                Vec3::new(boost * h, 0., 0.),
                Vec3::new(boost, 0., 0.),
            );
            state_near(out.step.states_after[i], expected, 3e-12);
        }
        near(out.actuator_work_j, baseline.actuator_work_j, 3e-13);
        near(
            out.constraint_work_defect_j,
            baseline.constraint_work_defect_j,
            3e-13,
        );
        let mut independently_paired = 0.;
        for pair in out
            .step
            .position_impulses
            .iter()
            .chain(&out.step.velocity_impulses)
        {
            let a = pair.first;
            let b = pair.second;
            assert_ne!(
                b.state_after.linear_momentum_world().x,
                b.state_before.linear_momentum_world().x
            );
            let before = b
                .point_kinematics_before
                .center_of_mass_velocity_world
                .sub(a.point_kinematics_before.center_of_mass_velocity_world);
            let after = b
                .point_kinematics_after
                .center_of_mass_velocity_world
                .sub(a.point_kinematics_after.center_of_mass_velocity_world);
            independently_paired += b.impulse_world.dot(before.scale(0.5).add(after.scale(0.5)));
            for event in [a, b] {
                independently_paired += event.angular_impulse_body.dot(
                    event
                        .point_kinematics_before
                        .angular_velocity_body
                        .scale(0.5)
                        .add(
                            event
                                .point_kinematics_after
                                .angular_velocity_body
                                .scale(0.5),
                        ),
                );
            }
        }
        near(
            out.actuator_work_j + out.constraint_work_defect_j,
            independently_paired,
            3e-15,
        );
    }
}

#[test]
fn g1_held_external_forces_keep_actuator_work_visible_under_large_boosts() {
    let p = [spherical(1.), spherical(2.), spherical(3.)];
    let links = [center_link(0, 1, 1.), center_link(1, 2, 2.)];
    let rates = [0.3, -0.2];
    let accelerations = [0.4, 0.1];
    let velocities = chain_rates(rates);
    let forces = [0.7, -0.2, 0.1];
    let acceleration = chain_rates(accelerations).map(|a| a + 0.1);
    let commands = std::array::from_fn::<_, 2, _>(|i| DistanceDrive {
        range_rate_mps: rates[i],
        range_acceleration_mps2: accelerations[i],
    });
    let wrenches = forces.map(|force| Wrench {
        force_world: Vec3::new(force, 0., 0.),
        torque_body: Vec3::ZERO,
    });
    let h = 0.2;
    let impulses = [
        0.5 * h * (forces[0] - acceleration[0]),
        0.5 * h * (3. * acceleration[2] - forces[2]),
    ];
    let expected_source: f64 = (0..2)
        .map(|i| impulses[i] * (2. * rates[i] + h * accelerations[i]))
        .sum();
    let mut values = Vec::new();
    for boost in [0., 2_f64.powi(30)] {
        let s = std::array::from_fn::<_, 3, _>(|i| {
            state(
                p[i],
                Vec3::new([0., 1., 3.][i], 0., 0.),
                UnitQuaternion::IDENTITY,
                Vec3::new(velocities[i], boost, 0.),
                Vec3::ZERO,
            )
        });
        let out = advance_driven_distance_graph(
            &s,
            &p,
            &links,
            &commands,
            free(),
            &wrenches,
            h,
            controls(),
            || false,
        )
        .unwrap();
        near(out.actuator_work_j, expected_source, 3e-14);
        near(out.constraint_work_defect_j, 0., 3e-14);
        let mut kinetic_change = 0.;
        let mut applied_work = 0.;
        for i in 0..3 {
            let next = out.step.states_after[i];
            let v = next.center_of_mass_velocity_world(p[i]).unwrap();
            near(v.x, velocities[i] + h * acceleration[i], 3e-14);
            near(v.y, boost, 0.);
            kinetic_change += 0.5 * p[i].mass() * (v.x - velocities[i]) * (v.x + velocities[i]);
            applied_work +=
                forces[i] * (next.pose().position_world().x - s[i].pose().position_world().x);
        }
        near(
            kinetic_change,
            applied_work + out.actuator_work_j + out.constraint_work_defect_j,
            3e-14,
        );
        if boost != 0. {
            // Subtracting the enormous laboratory-frame energies would lose
            // this physical work; the native relative account must not do so.
            near(
                totals(&out.step.states_after, &p).2 - totals(&s, &p).2,
                0.,
                0.,
            );
            assert!(out.actuator_work_j > 0.01);
        }
        values.push((out.actuator_work_j, out.constraint_work_defect_j));
    }
    assert_eq!(values[0], values[1]);
}

fn off_axis() -> ([MassProperties; 3], [RigidBodyState; 3], [DistanceLink; 3]) {
    let p = std::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1. + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1., 1.2).scale(1. + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.3 + 0.2 * i as f64).unwrap(),
        )
        .unwrap()
    });
    let positions = [
        Vec3::new(-0.8, 0.2, 0.1),
        Vec3::new(0.6, 0.1, -0.2),
        Vec3::new(0.3, 1.3, 0.4),
    ];
    let omega = Vec3::new(0.2, -0.1, 0.3);
    let s = std::array::from_fn(|i| {
        let q =
            UnitQuaternion::from_axis_angle(Vec3::new(1., -0.5, 2.), 0.2 - 0.3 * i as f64).unwrap();
        state(
            p[i],
            positions[i],
            q,
            Vec3::new(0.1, -0.2, 0.3).add(omega.cross(positions[i])),
            q.rotate_world_to_body(omega),
        )
    });
    let arms = [
        [Vec3::new(0.2, -0.1, 0.1), Vec3::new(-0.1, 0.2, -0.05)],
        [Vec3::new(0.1, 0.05, -0.15), Vec3::new(-0.2, -0.1, 0.1)],
        [Vec3::new(0.05, 0.15, 0.1), Vec3::new(0.1, -0.1, 0.05)],
    ];
    let links = std::array::from_fn(|i| {
        let bodies = [[0, 1], [1, 2], [2, 0]][i];
        let a = s[bodies[0]]
            .pose()
            .point_world_from_body(arms[i][0])
            .unwrap();
        let b = s[bodies[1]]
            .pose()
            .point_world_from_body(arms[i][1])
            .unwrap();
        DistanceLink {
            bodies,
            points_body_m: arms[i],
            length_m: b.sub(a).norm_squared().sqrt(),
        }
    });
    (p, s, links)
}
fn acceleration_commands() -> [DistanceDrive; 3] {
    [0.1, -0.2, 0.15].map(|a| DistanceDrive {
        range_rate_mps: 0.,
        range_acceleration_mps2: a,
    })
}

#[test]
fn g0_full_tensor_closed_actuated_framework_retains_momenta_work_and_reversibility() {
    let (p, initial, initial_links) = off_axis();
    let before = totals(&initial, &p);
    let mut s = initial.to_vec();
    let mut links = initial_links;
    let mut drives = acceleration_commands();
    let mut work = 0.;
    let mut source = 0.;
    let mut defect = 0.;
    for _ in 0..100 {
        let out = run(&s, &p, &links, &drives, 0.005);
        for event in out
            .step
            .position_impulses
            .iter()
            .chain(&out.step.velocity_impulses)
        {
            vector(event.impulse_balance_world, Vec3::ZERO, 0.);
            vector(event.angular_impulse_world, Vec3::ZERO, 5e-15);
        }
        for i in 0..3 {
            links[i] = out.links[i].link_after;
            drives[i] = out.links[i].drive_after;
            assert_eq!(
                out.step.states_after[i].pose(),
                out.step.free_steps[i].state_after.pose()
            );
        }
        work += out.step.constraint_work_j;
        source += out.actuator_work_j;
        defect += out.constraint_work_defect_j;
        s = out.step.states_after;
    }
    let after = totals(&s, &p);
    vector(before.0, after.0, 2e-12);
    vector(before.1, after.1, 2e-11);
    near(after.2 - before.2, work, 3e-12);
    near(work, source + defect, 3e-12);
    assert!(
        defect.abs() > 1e-11,
        "finite numerical work remains visible: {defect:e}"
    );
    s = s.into_iter().map(reverse).collect();
    for d in &mut drives {
        d.range_rate_mps = -d.range_rate_mps;
    }
    for _ in 0..100 {
        let out = run(&s, &p, &links, &drives, 0.005);
        for i in 0..3 {
            links[i] = out.links[i].link_after;
            drives[i] = out.links[i].drive_after;
        }
        s = out.step.states_after;
    }
    for i in 0..3 {
        state_near(reverse(s[i]), initial[i], 3e-10);
        near(links[i].length_m, initial_links[i].length_m, 3e-13);
        near(drives[i].range_rate_mps, 0., 3e-15);
    }
}

#[test]
fn g1_rotating_telescopic_pair_matches_continuous_phase_and_actuator_work() {
    let p = [spherical(1.), spherical(2.)];
    let mu = 2. / 3.;
    let rate = 0.2;
    let spin = 0.5;
    let horizon = 1.;
    let radius = 1. + rate * horizon;
    let angle = spin * horizon / radius;
    let exact = Vec3::new(radius * f64::cos(angle), radius * f64::sin(angle), 0.);
    let exact_work = 0.5 * mu * spin * spin * (1. / (radius * radius) - 1.);
    let results = [25, 50, 100].map(|count| {
        let mut s = [-2. / 3., 1. / 3.].map(|f| {
            state(
                p[if f < 0. { 0 } else { 1 }],
                Vec3::new(f, 0., 0.),
                UnitQuaternion::IDENTITY,
                Vec3::new(rate, spin, 0.).scale(f),
                Vec3::ZERO,
            )
        });
        let initial = totals(&s, &p);
        let mut link = center_link(0, 1, 1.);
        let mut drive = DistanceDrive {
            range_rate_mps: rate,
            range_acceleration_mps2: 0.,
        };
        let mut source = 0.;
        let mut defect = 0.;
        for _ in 0..count {
            let out = run(&s, &p, &[link], &[drive], horizon / count as f64);
            s.copy_from_slice(&out.step.states_after);
            link = out.links[0].link_after;
            drive = out.links[0].drive_after;
            source += out.actuator_work_j;
            defect += out.constraint_work_defect_j;
        }
        let final_totals = totals(&s, &p);
        vector(final_totals.0, initial.0, 2e-14);
        vector(final_totals.1, initial.1, 2e-13);
        near(final_totals.2 - initial.2, exact_work, 2e-13);
        near(source + defect, exact_work, 2e-13);
        assert!(source < 0.);
        let delta = s[1]
            .pose()
            .position_world()
            .sub(s[0].pose().position_world())
            .sub(exact);
        (
            delta.norm_squared().sqrt(),
            (source - exact_work).abs(),
            defect.abs(),
        )
    });
    for window in results.windows(2) {
        for ratio in [
            window[0].0 / window[1].0,
            window[0].1 / window[1].1,
            window[0].2 / window[1].2,
        ] {
            assert!((3.8..4.2).contains(&ratio), "{results:?}, ratio {ratio}");
        }
    }
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [v, i, j, k] = b.components();
    UnitQuaternion::new(
        w * v - x * i - y * j - z * k,
        w * i + x * v + y * k - z * j,
        w * j - x * k + y * v + z * i,
        w * k + x * j - y * i + z * v,
    )
    .unwrap()
}
fn transform(
    s: RigidBodyState,
    p: MassProperties,
    q: UnitQuaternion,
    shift: Vec3,
    boost: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            q.rotate_body_to_world(s.pose().position_world()).add(shift),
            compose(q, s.pose().orientation()),
        )
        .unwrap(),
        q.rotate_body_to_world(s.linear_momentum_world())
            .add(boost.scale(p.mass())),
        s.angular_momentum_body(),
    )
    .unwrap()
}

#[test]
fn g3_actuation_is_frame_boost_and_graph_order_covariant() {
    let (p, s, links) = off_axis();
    let drives = acceleration_commands();
    let h = 0.03;
    let base = run(&s, &p, &links, &drives, h);
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.8).unwrap();
    let shift = Vec3::new(3., -2., 4.);
    let boost = Vec3::new(0.7, -0.4, 0.2);
    let transformed = std::array::from_fn::<_, 3, _>(|i| transform(s[i], p[i], q, shift, boost));
    let out = run(&transformed, &p, &links, &drives, h);
    for i in 0..3 {
        state_near(
            out.step.states_after[i],
            transform(
                base.step.states_after[i],
                p[i],
                q,
                shift.add(boost.scale(h)),
                boost,
            ),
            2e-11,
        );
    }
    near(out.actuator_work_j, base.actuator_work_j, 3e-13);
    near(
        out.constraint_work_defect_j,
        base.constraint_work_defect_j,
        3e-13,
    );
    for order in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let rods = order.map(|i| DistanceLink {
            bodies: [links[i].bodies[1], links[i].bodies[0]],
            points_body_m: [links[i].points_body_m[1], links[i].points_body_m[0]],
            ..links[i]
        });
        let commands = order.map(|i| drives[i]);
        let out = run(&s, &p, &rods, &commands, h);
        for i in 0..3 {
            state_near(out.step.states_after[i], base.step.states_after[i], 2e-11);
        }
        near(out.actuator_work_j, base.actuator_work_j, 3e-13);
        near(
            out.constraint_work_defect_j,
            base.constraint_work_defect_j,
            3e-13,
        );
    }
}

#[test]
fn g0_polynomial_admission_and_real_maximum_graph_honor_resource_bounds() {
    let p = vec![spherical(1.); 33];
    let s: Vec<_> = (0..33)
        .map(|i| {
            state(
                p[i],
                Vec3::new(i as f64, 0., 0.),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
                Vec3::ZERO,
            )
        })
        .collect();
    let links: Vec<_> = (0..32).map(|i| center_link(i, i + 1, 1.)).collect();
    let drives: Vec<_> = (0..32)
        .map(|i| DistanceDrive {
            range_rate_mps: 0.,
            range_acceleration_mps2: if i % 2 == 0 { 0.01 } else { -0.01 },
        })
        .collect();
    let out = run(&s, &p, &links, &drives, 0.1);
    assert_eq!(out.links.len(), 32);
    assert!(out.actuator_work_j > 0.);
    let mut cap = controls();
    cap.maximum_live_bytes = out.planned_live_bytes;
    let admitted = advance_driven_distance_graph(
        &s,
        &p,
        &links,
        &drives,
        free(),
        &vec![Wrench::ZERO; 33],
        0.1,
        cap,
        || false,
    )
    .unwrap();
    assert_eq!(admitted, out);
    cap.maximum_live_bytes -= 1;
    assert!(matches!(
        advance_driven_distance_graph(
            &s,
            &p,
            &links,
            &drives,
            free(),
            &vec![Wrench::ZERO; 33],
            0.1,
            cap,
            || false
        ),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    let more = vec![center_link(0, 1, 1.); 33];
    assert!(matches!(
        advance_driven_distance_graph(
            &s,
            &p,
            &more,
            &vec![DistanceDrive::default(); 33],
            free(),
            &vec![Wrench::ZERO; 33],
            0.1,
            controls(),
            || false
        ),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));
    // Positive endpoint lengths do not admit a path passing through zero.
    for (length, rate, acceleration, h) in [
        (1., -4., 8., 1.),
        (1., -2., 0., 1.),
        (2_f64.powi(20), 2_f64.powi(-40), 0., 1.),
        (2_f64.powi(22), 2_f64.powi(20), 2_f64.powi(-34), 1.),
        (1., 1e308, 1e308, 1e308),
    ] {
        let rod = center_link(0, 1, length);
        let drive = DistanceDrive {
            range_rate_mps: rate,
            range_acceleration_mps2: acceleration,
        };
        assert!(
            matches!(
                advance_driven_distance_graph(
                    &s[..2],
                    &p[..2],
                    &[rod],
                    &[drive],
                    free(),
                    &[Wrench::ZERO; 2],
                    h,
                    controls(),
                    || false
                ),
                Err(DistanceGraphError::PrescribedMotion { constraint: 0, .. })
            ),
            "unadmitted prescription {length} {rate} {acceleration} {h}"
        );
    }
    assert!(matches!(
        advance_driven_distance_graph(
            &s[..2],
            &p[..2],
            &links[..1],
            &[DistanceDrive {
                range_rate_mps: 0.1,
                range_acceleration_mps2: 0.
            }],
            free(),
            &[Wrench::ZERO; 2],
            0.1,
            controls(),
            || false
        ),
        Err(DistanceGraphError::InitialConstraint { constraint: 0 })
    ));
}

#[test]
fn g4_every_cancellation_boundary_preserves_inputs_and_retries_the_same_step() {
    let (p, s, links) = off_axis();
    let drives = acceleration_commands();
    let mut calls = 0;
    let reference = advance_driven_distance_graph(
        &s,
        &p,
        &links,
        &drives,
        free(),
        &[Wrench::ZERO; 3],
        0.02,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in 1..=calls {
        let mut n = 0;
        let result = advance_driven_distance_graph(
            &s,
            &p,
            &links,
            &drives,
            free(),
            &[Wrench::ZERO; 3],
            0.02,
            controls(),
            || {
                n += 1;
                n == stop
            },
        );
        assert!(
            matches!(
                result,
                Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
            ),
            "checkpoint {stop}/{calls}: {result:?}"
        );
        assert_eq!(s, off_axis().1);
    }
    assert_eq!(run(&s, &p, &links, &drives, 0.02), reference);
}
