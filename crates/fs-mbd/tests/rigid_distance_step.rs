use fs_mbd::rigid_rattle::distance::{
    DistanceStepControls, DistanceStepError, RigidDistanceConstraint,
};
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
fn controls() -> DistanceStepControls {
    DistanceStepControls {
        maximum_impulse_ns: 10.0,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 40,
        maximum_backtracks: 16,
    }
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap())
}
fn spherical(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}
fn state(
    properties: MassProperties,
    position: Vec3,
    orientation: UnitQuaternion,
    velocity: Vec3,
    omega: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn constraint(
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    length: f64,
) -> RigidDistanceConstraint {
    RigidDistanceConstraint::new(properties, arms, length, 2e-13, 2e-13).unwrap()
}
fn totals(states: [RigidBodyState; 2], properties: [MassProperties; 2]) -> (Vec3, Vec3, f64) {
    let d = [
        free().diagnostics(states[0], properties[0]).unwrap(),
        free().diagnostics(states[1], properties[1]).unwrap(),
    ];
    (
        d[0].linear_momentum_world.add(d[1].linear_momentum_world),
        d[0].angular_momentum_world.add(d[1].angular_momentum_world),
        d[0].mechanical_energy + d[1].mechanical_energy,
    )
}
fn dimer(speed: f64, boost: Vec3) -> ([MassProperties; 2], [RigidBodyState; 2]) {
    let p = [spherical(2.0), spherical(3.0)];
    let s = [
        state(
            p[0],
            Vec3::new(-1.2, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            boost.add(Vec3::new(0.0, -0.6 * speed, 0.0)),
            Vec3::ZERO,
        ),
        state(
            p[1],
            Vec3::new(0.8, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            boost.add(Vec3::new(0.0, 0.4 * speed, 0.0)),
            Vec3::ZERO,
        ),
    ];
    (p, s)
}

#[test]
fn unequal_mass_dimer_matches_closed_form_pair_recoil_and_orbit() {
    let speed = 1.5;
    let boost = Vec3::new(0.3, -0.2, 0.4);
    let (p, s) = dimer(speed, boost);
    let dt = 0.1;
    let out = constraint(p, [Vec3::ZERO; 2], 2.0)
        .advance(s, free(), [Wrench::ZERO; 2], dt, controls(), || false)
        .unwrap();
    // Local root of |(L+h*W*lambda,h*v)|=L, with W=1/m0+1/m1.
    let x = (4.0 - dt * dt * speed * speed).sqrt();
    let y = dt * speed;
    let mobility = 1.0 / 2.0 + 1.0 / 3.0;
    let expected_impulse = -dt * speed * speed / (mobility * (2.0 + x));
    vector(
        out.position_impulse.second.impulse_world,
        Vec3::new(expected_impulse, 0.0, 0.0),
        2e-12,
    );
    let relative = Vec3::new(x, y, 0.0);
    let relative_velocity = Vec3::new(-speed * y / 2.0, speed * x / 2.0, 0.0);
    for (i, factor) in [-0.6, 0.4].into_iter().enumerate() {
        vector(
            out.states_after[i].pose().position_world(),
            boost.scale(dt).add(relative.scale(factor)),
            2e-13,
        );
        vector(
            out.states_after[i]
                .center_of_mass_velocity_world(p[i])
                .unwrap(),
            boost.add(relative_velocity.scale(factor)),
            2e-12,
        );
        assert_eq!(
            out.states_after[i].pose(),
            out.free_steps[i].state_after.pose()
        );
    }
    let before = totals(s, p);
    let after = totals(out.states_after, p);
    vector(after.0, before.0, 2e-15);
    vector(after.1, before.1, 3e-13);
    near(after.2, before.2, 2e-12);
    near(out.final_gap_m, 0.0, 2e-13);
    near(out.final_range_rate_mps, 0.0, 2e-13);
}

#[test]
fn unequal_external_loads_share_acceleration_without_pinning_either_body() {
    let (p, s) = dimer(0.0, Vec3::ZERO);
    let w = [
        Wrench {
            force_world: Vec3::new(1.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(4.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
    ];
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let stepper = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap());
    let dt = 0.08;
    let out = constraint(p, [Vec3::ZERO; 2], 2.0)
        .advance(s, stepper, w, dt, controls(), || false)
        .unwrap();
    let acceleration = Vec3::new(1.0, 0.0, -9.81);
    for i in 0..2 {
        vector(
            out.states_after[i].pose().position_world(),
            s[i].pose()
                .position_world()
                .add(acceleration.scale(0.5 * dt * dt)),
            2e-14,
        );
        vector(
            out.states_after[i]
                .center_of_mass_velocity_world(p[i])
                .unwrap(),
            acceleration.scale(dt),
            2e-14,
        );
    }
    let expected = -0.5 * dt * ((4.0 / 3.0) - (1.0 / 2.0)) / (1.0 / 2.0 + 1.0 / 3.0);
    near(out.position_impulse.second.impulse_world.x, expected, 2e-13);
    near(out.velocity_impulse.second.impulse_world.x, expected, 2e-13);
    vector(
        totals(out.states_after, p).0,
        acceleration.scale(5.0 * dt),
        3e-14,
    );
}

fn off_axis() -> ([MassProperties; 2], [Vec3; 2], [RigidBodyState; 2]) {
    let p = [
        MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.7, 1.0, 1.2),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3.0,
            Vec3::ZERO,
            Vec3::new(0.8, 1.3, 1.4),
            UnitQuaternion::from_axis_angle(Vec3::new(-1.0, 1.0, 2.0), -0.7).unwrap(),
        )
        .unwrap(),
    ];
    let q = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 0.4), 0.4).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-0.5, 1.0, 1.0), -0.5).unwrap(),
    ];
    let arms = [Vec3::new(0.3, -0.15, 0.2), Vec3::new(-0.2, 0.25, 0.15)];
    let world = [
        q[0].rotate_body_to_world(arms[0]),
        q[1].rotate_body_to_world(arms[1]),
    ];
    let omega = [Vec3::new(0.3, -0.2, 0.4), Vec3::new(-0.2, 0.1, 0.3)];
    let spin = [
        q[0].rotate_body_to_world(omega[0]).cross(world[0]),
        q[1].rotate_body_to_world(omega[1]).cross(world[1]),
    ];
    let a = Vec3::new(-0.6, 0.2, -0.1);
    let b = a.add(world[0]).add(Vec3::new(1.5, 0.0, 0.0)).sub(world[1]);
    let va = Vec3::new(0.1, 0.2, -0.1);
    let vb = va.add(spin[0]).sub(spin[1]).add(Vec3::new(0.0, 0.4, 0.2));
    let s = [
        state(p[0], a, q[0], va, omega[0]),
        state(p[1], b, q[1], vb, omega[1]),
    ];
    (p, arms, s)
}

#[test]
fn full_inertia_offset_points_conserve_all_linear_and_angular_momentum() {
    let (p, arms, mut s) = off_axis();
    let law = constraint(p, arms, 1.5);
    let initial = totals(s, p);
    let mut work = 0.0;
    for _ in 0..300 {
        let out = law
            .advance(s, free(), [Wrench::ZERO; 2], 0.005, controls(), || false)
            .unwrap();
        for event in [&out.position_impulse, &out.velocity_impulse] {
            vector(event.impulse_balance_world, Vec3::ZERO, 0.0);
            vector(event.angular_impulse_world, Vec3::ZERO, 2e-15);
            near(event.work_energy_residual_j, 0.0, 2e-16);
        }
        for i in 0..2 {
            assert_eq!(
                out.states_after[i].pose(),
                out.free_steps[i].state_after.pose()
            );
        }
        s = out.states_after;
        work += out.constraint_work_j;
        near(out.final_gap_m, 0.0, 2e-13);
        near(out.final_range_rate_mps, 0.0, 2e-13);
    }
    let final_totals = totals(s, p);
    vector(final_totals.0, initial.0, 2e-12);
    vector(final_totals.1, initial.1, 3e-12);
    near(final_totals.2 - initial.2, work, 2e-12);
    assert!(
        work.abs() > 1e-10,
        "nonzero finite-step work must remain visible"
    );
}

fn reversed(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.0),
        s.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}
fn state_near(a: RigidBodyState, b: RigidBodyState, tol: f64) {
    vector(a.pose().position_world(), b.pose().position_world(), tol);
    vector(a.linear_momentum_world(), b.linear_momentum_world(), tol);
    vector(a.angular_momentum_body(), b.angular_momentum_body(), tol);
    for e in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        vector(
            a.pose().orientation().rotate_body_to_world(e),
            b.pose().orientation().rotate_body_to_world(e),
            tol,
        );
    }
}

#[test]
fn local_branch_reverses_and_has_second_order_orbital_phase() {
    let (p, arms, initial) = off_axis();
    let law = constraint(p, arms, 1.5);
    let mut s = initial;
    for _ in 0..60 {
        s = law
            .advance(s, free(), [Wrench::ZERO; 2], 0.01, controls(), || false)
            .unwrap()
            .states_after;
    }
    s = s.map(reversed);
    for _ in 0..60 {
        s = law
            .advance(s, free(), [Wrench::ZERO; 2], 0.01, controls(), || false)
            .unwrap()
            .states_after;
    }
    for i in 0..2 {
        state_near(reversed(s[i]), initial[i], 2e-10);
    }

    let (p, initial) = dimer(1.4, Vec3::ZERO);
    let law = constraint(p, [Vec3::ZERO; 2], 2.0);
    let errors = [10, 20, 40].map(|steps| {
        let mut s = initial;
        for _ in 0..steps {
            s = law
                .advance(
                    s,
                    free(),
                    [Wrench::ZERO; 2],
                    1.0 / f64::from(steps),
                    controls(),
                    || false,
                )
                .unwrap()
                .states_after;
        }
        let exact = Vec3::new(0.8 * 0.7_f64.cos(), 0.8 * 0.7_f64.sin(), 0.0);
        s[1].pose()
            .position_world()
            .sub(exact)
            .norm_squared()
            .sqrt()
    });
    for pair in errors.windows(2) {
        assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
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
    r: UnitQuaternion,
    shift: Vec3,
    boost: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            r.rotate_body_to_world(s.pose().position_world()).add(shift),
            compose(r, s.pose().orientation()),
        )
        .unwrap(),
        r.rotate_body_to_world(s.linear_momentum_world())
            .add(boost.scale(p.mass())),
        s.angular_momentum_body(),
    )
    .unwrap()
}

#[test]
fn pair_motion_is_objective_boost_covariant_and_body_order_independent() {
    let (p, arms, initial) = off_axis();
    let law = constraint(p, arms, 1.5);
    let dt = 0.03;
    let out = law
        .advance(initial, free(), [Wrench::ZERO; 2], dt, controls(), || false)
        .unwrap();
    let r = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.8).unwrap();
    let shift = Vec3::new(2.0, -4.0, 1.0);
    let boost = Vec3::new(0.7, -0.4, 0.3);
    let moved = std::array::from_fn(|i| transform(initial[i], p[i], r, shift, boost));
    let output = law
        .advance(moved, free(), [Wrench::ZERO; 2], dt, controls(), || false)
        .unwrap();
    for i in 0..2 {
        let expected = transform(
            out.states_after[i],
            p[i],
            r,
            shift.add(boost.scale(dt)),
            boost,
        );
        state_near(output.states_after[i], expected, 2e-12);
    }
    let swapped = constraint([p[1], p[0]], [arms[1], arms[0]], 1.5)
        .advance(
            [initial[1], initial[0]],
            free(),
            [Wrench::ZERO; 2],
            dt,
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..2 {
        state_near(out.states_after[i], swapped.states_after[1 - i], 2e-14);
    }
}

#[test]
fn geometry_local_limits_and_every_cancellation_boundary_refuse_atomically() {
    let (p, s) = dimer(1.5, Vec3::ZERO);
    let law = constraint(p, [Vec3::ZERO; 2], 2.0);
    assert!(RigidDistanceConstraint::new(p, [Vec3::ZERO; 2], 0.0, 1e-12, 1e-12).is_err());
    assert!(RigidDistanceConstraint::new(p, [Vec3::ZERO; 2], 1.0, 1.0, 1e-12).is_err());
    let collapsed = [
        s[0],
        RigidBodyState::new(s[0].pose(), s[1].linear_momentum_world(), Vec3::ZERO).unwrap(),
    ];
    assert!(matches!(
        law.advance(
            collapsed,
            free(),
            [Wrench::ZERO; 2],
            0.1,
            controls(),
            || false
        ),
        Err(DistanceStepError::Invalid(_))
    ));
    let mut fast = s;
    fast[0] = RigidBodyState::new(
        s[0].pose(),
        s[0].linear_momentum_world().add(Vec3::new(0.1, 0.0, 0.0)),
        Vec3::ZERO,
    )
    .unwrap();
    assert_eq!(
        law.advance(fast, free(), [Wrench::ZERO; 2], 0.1, controls(), || false),
        Err(DistanceStepError::InitialConstraint)
    );
    let mut cap = controls();
    cap.maximum_rotation_increment_rad = 0.01;
    assert_eq!(
        law.advance(s, free(), [Wrench::ZERO; 2], 0.1, cap, || false),
        Err(DistanceStepError::LocalBranch)
    );
    cap = controls();
    cap.maximum_impulse_ns = 1e-12;
    assert!(
        law.advance(s, free(), [Wrench::ZERO; 2], 0.1, cap, || false)
            .is_err()
    );
    cap = controls();
    cap.maximum_position_iterations = 1;
    assert_eq!(
        law.advance(s, free(), [Wrench::ZERO; 2], 0.1, cap, || false),
        Err(DistanceStepError::NotConverged)
    );
    let mut polls = 0;
    let accepted = law
        .advance(s, free(), [Wrench::ZERO; 2], 0.1, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert_eq!(
            law.advance(s, free(), [Wrench::ZERO; 2], 0.1, controls(), || {
                count += 1;
                count == stop
            }),
            Err(DistanceStepError::Cancelled)
        );
    }
    assert_eq!(
        accepted,
        law.advance(s, free(), [Wrench::ZERO; 2], 0.1, controls(), || false)
            .unwrap()
    );
    vector(s[0].pose().position_world(), Vec3::new(-1.2, 0.0, 0.0), 0.0);
}
