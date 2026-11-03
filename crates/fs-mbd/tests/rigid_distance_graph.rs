use fs_mbd::rigid_rattle::distance::graph::{
    DistanceGraphControls, DistanceGraphError, DistanceGraphStep, DistanceLink,
    advance_distance_graph,
};
use fs_mbd::rigid_rattle::distance::{
    DistanceStepControls, DistanceStepError, RigidDistanceConstraint,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn controls() -> DistanceGraphControls {
    DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 10.0,
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
fn resting(p: MassProperties, x: Vec3) -> RigidBodyState {
    state(p, x, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO)
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
    rods: &[DistanceLink],
    h: f64,
) -> DistanceGraphStep {
    advance_distance_graph(
        s,
        p,
        rods,
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
        (Vec3::ZERO, Vec3::ZERO, 0.0),
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
        s.linear_momentum_world().scale(-1.0),
        s.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}
fn chain() -> (
    [MassProperties; 3],
    [RigidBodyState; 3],
    [DistanceLink; 2],
    [Wrench; 3],
) {
    let p = [spherical(1.), spherical(2.), spherical(3.)];
    let s = std::array::from_fn(|i| resting(p[i], Vec3::new([0., 1., 3.][i], 0., 0.)));
    let w = [7., -2., 1.].map(|f| Wrench {
        force_world: Vec3::new(f, 0., 0.),
        torque_body: Vec3::ZERO,
    });
    (p, s, [center_link(0, 1, 1.), center_link(1, 2, 2.)], w)
}

#[test]
fn g1_shared_body_chain_uses_all_cross_responses_and_physical_recoil() {
    let (p, s, rods, w) = chain();
    let h = 0.1;
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let stepper = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap());
    let out = advance_distance_graph(&s, &p, &rods, stepper, &w, h, controls(), || false).unwrap();
    let acceleration = Vec3::new(1.0, 0.0, -9.81);
    for i in 0..3 {
        vector(
            out.states_after[i].pose().position_world(),
            s[i].pose()
                .position_world()
                .add(acceleration.scale(0.5 * h * h)),
            3e-14,
        );
        vector(
            out.states_after[i]
                .center_of_mass_velocity_world(p[i])
                .unwrap(),
            acceleration.scale(h),
            3e-14,
        );
        assert_eq!(
            out.states_after[i].pose(),
            out.free_steps[i].state_after.pose()
        );
    }
    // From per-body momenta: -j0 = h(m0*a-F0)/2 and j1=h(m2*a-F2)/2.
    for events in [&out.position_impulses, &out.velocity_impulses] {
        near(events[0].second.impulse_world.x, 0.3, 3e-14);
        near(events[1].second.impulse_world.x, 0.1, 3e-14);
    }
    // Individual pair work depends on application order on the middle body;
    // the sum must include those cross terms and equal the actual change.
    let event_work: f64 = out
        .position_impulses
        .iter()
        .chain(&out.velocity_impulses)
        .map(|e| e.impulse_work_j)
        .sum();
    near(out.constraint_work_j, event_work, 2e-16);
    let free_change: f64 = out
        .free_steps
        .iter()
        .map(|s| s.diagnostics_after.mechanical_energy - s.diagnostics_before.mechanical_energy)
        .sum();
    let before: f64 = s
        .iter()
        .zip(p)
        .map(|(&s, p)| stepper.diagnostics(s, p).unwrap().mechanical_energy)
        .sum();
    let after: f64 = out
        .states_after
        .iter()
        .zip(p)
        .map(|(&s, p)| stepper.diagnostics(s, p).unwrap().mechanical_energy)
        .sum();
    near(after - before, free_change + out.constraint_work_j, 2e-14);
    for r in out.residuals {
        near(r.gap_m, 0., 2e-13);
        near(r.range_rate_mps, 0., 2e-13);
    }
}

fn triangle() -> ([MassProperties; 3], [RigidBodyState; 3], [DistanceLink; 3]) {
    let p = [spherical(1.); 3];
    let omega = 0.7;
    let s = std::array::from_fn(|i| {
        let angle = std::f64::consts::TAU * i as f64 / 3.;
        let x = Vec3::new(angle.cos(), angle.sin(), 0.0);
        state(
            p[i],
            x,
            UnitQuaternion::IDENTITY,
            Vec3::new(0., 0., omega).cross(x),
            Vec3::ZERO,
        )
    });
    let rods = [[0, 1], [1, 2], [2, 0]].map(|[a, b]| {
        center_link(
            a,
            b,
            s[b].pose()
                .position_world()
                .sub(s[a].pose().position_world())
                .norm_squared()
                .sqrt(),
        )
    });
    (p, s, rods)
}

#[test]
fn g1_closed_triangle_has_analytic_collective_rotation_and_second_order_phase() {
    let (p, initial, rods) = triangle();
    let omega = 0.7;
    let h = 0.1;
    let out = run(&initial, &p, &rods, h);
    let angle = (h * omega).asin();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0., 0., 1.), angle).unwrap();
    for i in 0..3 {
        vector(
            out.states_after[i].pose().position_world(),
            rotation.rotate_body_to_world(initial[i].pose().position_world()),
            5e-13,
        );
        vector(
            out.states_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(initial[i].linear_momentum_world()),
            8e-13,
        );
    }
    let errors = [10, 20, 40].map(|n| {
        let mut s = initial.to_vec();
        for _ in 0..n {
            s = run(&s, &p, &rods, 1. / f64::from(n)).states_after;
        }
        let exact = Vec3::new(0.7_f64.cos(), 0.7_f64.sin(), 0.);
        s[0].pose()
            .position_world()
            .sub(exact)
            .norm_squared()
            .sqrt()
    });
    for w in errors.windows(2) {
        assert!((3.9..4.1).contains(&(w[0] / w[1])), "{errors:?}");
    }
    let a = totals(&initial, &p);
    let b = totals(&out.states_after, &p);
    vector(a.0, b.0, 8e-15);
    vector(a.1, b.1, 3e-13);
    near(a.2, b.2, 3e-13);
}

fn off_axis() -> ([MassProperties; 3], [RigidBodyState; 3], [DistanceLink; 3]) {
    let p = std::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1.0 + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1.0, 1.2).scale(1.0 + 0.2 * i as f64),
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
    let rods = std::array::from_fn(|i| {
        let [a, b] = [[0, 1], [1, 2], [2, 0]][i];
        let x = s[a].pose().point_world_from_body(arms[i][0]).unwrap();
        let y = s[b].pose().point_world_from_body(arms[i][1]).unwrap();
        DistanceLink {
            bodies: [a, b],
            points_body_m: arms[i],
            length_m: y.sub(x).norm_squared().sqrt(),
        }
    });
    (p, s, rods)
}

#[test]
fn g0_full_inertia_shared_attachments_preserve_momenta_and_keep_constraint_work() {
    let (p, initial, rods) = off_axis();
    let before = totals(&initial, &p);
    let mut s = initial.to_vec();
    let mut work = 0.;
    for _ in 0..300 {
        let out = run(&s, &p, &rods, 0.005);
        for event in out.position_impulses.iter().chain(&out.velocity_impulses) {
            vector(event.angular_impulse_world, Vec3::ZERO, 5e-15);
            vector(event.impulse_balance_world, Vec3::ZERO, 0.);
            near(event.work_energy_residual_j, 0., 3e-16);
        }
        for i in 0..3 {
            assert_eq!(
                out.states_after[i].pose(),
                out.free_steps[i].state_after.pose()
            );
        }
        work += out.constraint_work_j;
        s = out.states_after;
    }
    let after = totals(&s, &p);
    vector(after.0, before.0, 2e-12);
    vector(after.1, before.1, 2e-11);
    near(after.2 - before.2, work, 3e-12);
    assert!(
        work.abs() > 1e-10,
        "finite-step work should remain observable: {work:e}"
    );
    let mut s = initial.to_vec();
    for _ in 0..40 {
        s = run(&s, &p, &rods, 0.01).states_after;
    }
    s = s.into_iter().map(reverse).collect();
    for _ in 0..40 {
        s = run(&s, &p, &rods, 0.01).states_after;
    }
    for i in 0..3 {
        state_near(reverse(s[i]), initial[i], 3e-10);
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
fn g3_frame_boost_body_and_constraint_permutations_preserve_the_same_motion() {
    let (p, s, rods) = off_axis();
    let dt = 0.03;
    let reference = run(&s, &p, &rods, dt);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1., 2., -1.), 0.8).unwrap();
    let shift = Vec3::new(3., -2., 4.);
    let boost = Vec3::new(0.7, -0.4, 0.2);
    let moved = std::array::from_fn::<_, 3, _>(|i| transform(s[i], p[i], rotation, shift, boost));
    let result = run(&moved, &p, &rods, dt);
    for i in 0..3 {
        state_near(
            result.states_after[i],
            transform(
                reference.states_after[i],
                p[i],
                rotation,
                shift.add(boost.scale(dt)),
                boost,
            ),
            2e-11,
        );
    }
    for order in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let edges = order.map(|i| rods[i]);
        let result = run(&s, &p, &edges, dt);
        for i in 0..3 {
            state_near(result.states_after[i], reference.states_after[i], 2e-11);
        }
        near(result.constraint_work_j, reference.constraint_work_j, 3e-13);
    }
    let order = [2, 0, 1];
    let mut inverse = [0; 3];
    for i in 0..3 {
        inverse[order[i]] = i;
    }
    let p2 = order.map(|i| p[i]);
    let s2 = order.map(|i| s[i]);
    let edges = rods.map(|e| DistanceLink {
        bodies: [inverse[e.bodies[1]], inverse[e.bodies[0]]],
        points_body_m: [e.points_body_m[1], e.points_body_m[0]],
        ..e
    });
    let result = run(&s2, &p2, &edges, dt);
    for i in 0..3 {
        state_near(
            result.states_after[inverse[i]],
            reference.states_after[i],
            2e-11,
        );
    }
}

#[test]
fn g1_pair_specialization_disconnected_islands_and_empty_graph_keep_native_owners() {
    let (p, s, rods) = off_axis();
    let dt = 0.02;
    let graph = run(&s, &p, &rods[..1], dt);
    let e = rods[0];
    let pair =
        RigidDistanceConstraint::new([p[0], p[1]], e.points_body_m, e.length_m, 2e-13, 2e-13)
            .unwrap()
            .advance(
                [s[0], s[1]],
                free(),
                [Wrench::ZERO; 2],
                dt,
                controls().step,
                || false,
            )
            .unwrap();
    for i in 0..2 {
        state_near(graph.states_after[i], pair.states_after[i], 8e-13);
    }
    let third = free().step(s[2], p[2], Wrench::ZERO, dt).unwrap();
    assert_eq!(graph.states_after[2], third.state_after);
    let empty = run(&s, &p, &[], dt);
    for i in 0..3 {
        assert_eq!(
            empty.states_after[i],
            free()
                .step(s[i], p[i], Wrench::ZERO, dt)
                .unwrap()
                .state_after
        );
    }
    assert!(empty.residuals.is_empty());
    assert_eq!(empty.constraint_work_j, 0.0);
    // Two disjoint copies of the same rod must match separate native pairs.
    let states = [s[0], s[1], s[0], s[1]];
    let props = [p[0], p[1], p[0], p[1]];
    let edges = [
        e,
        DistanceLink {
            bodies: [2, 3],
            ..e
        },
    ];
    let islands = run(&states, &props, &edges, dt);
    for i in 0..4 {
        state_near(islands.states_after[i], pair.states_after[i % 2], 8e-13);
    }
}

#[test]
fn g4_resource_regular_branch_and_shape_refusals_are_atomic() {
    let (p, s, rods, w) = chain();
    let baseline =
        advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, controls(), || false).unwrap();
    let exact = DistanceGraphControls {
        maximum_trial_evaluations: baseline.trial_evaluations,
        maximum_live_bytes: baseline.planned_live_bytes,
        ..controls()
    };
    assert_eq!(
        advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, exact, || false).unwrap(),
        baseline
    );
    for c in [
        DistanceGraphControls {
            maximum_trial_evaluations: baseline.trial_evaluations - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_live_bytes: baseline.planned_live_bytes - 1,
            ..exact
        },
        DistanceGraphControls {
            maximum_bodies: 2,
            ..exact
        },
        DistanceGraphControls {
            maximum_constraints: 1,
            ..exact
        },
    ] {
        assert!(matches!(
            advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, c, || false),
            Err(DistanceGraphError::Limit { .. })
        ));
    }
    assert!(matches!(
        advance_distance_graph(
            &s,
            &p,
            &[rods[0], rods[0]],
            free(),
            &w,
            0.1,
            controls(),
            || false
        ),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    let redundant = [rods[0], rods[1], center_link(0, 2, 3.)];
    assert!(matches!(
        advance_distance_graph(&s, &p, &redundant, free(), &w, 0.1, controls(), || false),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    for bad in [
        DistanceLink {
            bodies: [0, 0],
            ..rods[0]
        },
        DistanceLink {
            bodies: [0, 3],
            ..rods[0]
        },
        DistanceLink {
            length_m: 0.0,
            ..rods[0]
        },
    ] {
        assert!(
            advance_distance_graph(&s, &p, &[bad], free(), &w, 0.1, controls(), || false).is_err()
        );
    }
    let c = DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 1e-4,
            ..controls().step
        },
        ..controls()
    };
    assert_eq!(
        advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, c, || false),
        Err(DistanceGraphError::Step(DistanceStepError::LocalBranch))
    );
    assert!(
        advance_distance_graph(&s, &p[..2], &rods, free(), &w, 0.1, controls(), || false).is_err()
    );
    // The admitted 32-row boundary is a real connected framework, not empty work.
    let props = vec![spherical(1.); 34];
    let states: Vec<_> = (0..34)
        .map(|i| resting(props[i], Vec3::new(i as f64, 0., 0.)))
        .collect();
    let edges: Vec<_> = (0..33).map(|i| center_link(i, i + 1, 1.)).collect();
    let at_bound = run(&states, &props, &edges[..32], 0.01);
    assert_eq!(at_bound.states_after, states);
    assert_eq!(at_bound.residuals.len(), 32);
    assert!(matches!(
        advance_distance_graph(
            &states,
            &props,
            &edges,
            free(),
            &vec![Wrench::ZERO; 34],
            0.01,
            controls(),
            || false
        ),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            ..
        })
    ));
    assert_eq!(s, chain().1);
}

#[test]
fn g4_every_observed_checkpoint_refuses_without_partial_publication_and_retries() {
    let (p, s, rods, w) = chain();
    let mut checkpoints = 0;
    let expected = advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, controls(), || {
        checkpoints += 1;
        false
    })
    .unwrap();
    assert!(checkpoints > 20);
    for stop in 1..=checkpoints {
        let mut calls = 0;
        let result = advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, controls(), || {
            calls += 1;
            calls == stop
        });
        assert_eq!(
            result,
            Err(DistanceGraphError::Step(DistanceStepError::Cancelled)),
            "stop={stop}"
        );
    }
    assert_eq!(s, chain().1);
    assert_eq!(
        advance_distance_graph(&s, &p, &rods, free(), &w, 0.1, controls(), || false).unwrap(),
        expected
    );
}
