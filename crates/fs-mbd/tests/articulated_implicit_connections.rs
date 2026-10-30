//! Actual finite-attachment dynamics: stiff modes, coupled trees and spatial work.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, MAX_IMPLICIT_SPATIAL_EVALUATIONS,
    step_free_floating_spatial_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad, step_free_floating_spatial,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

const ZERO: Wrench = Wrench {
    torque: Vec3 {
        x: 0.,
        y: 0.,
        z: 0.,
    },
    force: Vec3 {
        x: 0.,
        y: 0.,
        z: 0.,
    },
};

fn close(a: f64, b: f64, eps: f64) {
    assert!(
        (a - b).abs() <= eps,
        "{a:.16e} != {b:.16e}, error {:.4e}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, eps: f64) {
    assert!(
        (a - b).norm() <= eps,
        "{a:?} != {b:?}, error {:.4e}",
        (a - b).norm()
    );
}
fn rot(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn pose(v: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), v).unwrap()
}
fn inertia(m: f64, com: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        m,
        com,
        Mat3 {
            m: [2., 0.2, 0.1, 0.2, 2.5, 0.15, 0.1, 0.15, 3.],
        },
    )
    .unwrap()
}
fn point(link: usize, p: Vec3) -> MaterialPoint {
    MaterialPoint {
        link,
        point_body: p,
    }
}
fn connection(a: usize, b: usize, k: f64, c: f64) -> AxialConnection {
    AxialConnection {
        a: point(a, Vec3::default()),
        b: point(b, Vec3::default()),
        rest_length_m: 1.,
        stiffness_n_per_m: k,
        damping_ns_per_m: c,
    }
}
fn state(q: Vec<f64>, v: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: q,
        joint_velocity: v,
    }
}
fn balance(step: &FreeFloatingMotionStep, eps: f64) {
    vector(
        step.momentum_after.linear_momentum_world - step.momentum_before.linear_momentum_world,
        step.external_linear_impulse_world,
        eps,
    );
    vector(
        step.momentum_after.angular_momentum_world - step.momentum_before.angular_momentum_world,
        step.external_angular_impulse_world,
        eps,
    );
}
fn slider() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3., Vec3::default()),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(Vec3::new(1., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2., Vec3::default()),
        ),
    ])
    .unwrap()
}
fn held<'a>(effort: &'a [f64], wrench: &'a [Wrench]) -> FreeFloatingMotionLoads<'a> {
    FreeFloatingMotionLoads {
        generalized_force: effort,
        external_body_wrench: wrench,
        gravity_world: Vec3::default(),
    }
}
fn midpoint_mode(
    q: f64,
    v: f64,
    stiffness_per_mass: f64,
    damping_per_mass: f64,
    h: f64,
) -> (f64, f64, f64) {
    // Independently solve the two scalar trapezoidal equations, without native ABA.
    let mid = (v - 0.5 * h * stiffness_per_mass * q)
        / (1. + 0.5 * h * damping_per_mass + 0.25 * h * h * stiffness_per_mass);
    (q + h * mid, 2. * mid - v, mid)
}

#[test]
fn unresolved_stiff_linear_modes_remain_bounded_and_keep_physical_recoil_and_loss() {
    let model = slider();
    let mu = 1.2;
    let h = 0.01;
    for damping in [0., 240.] {
        let springs = [connection(0, 1, 1.2e6, damping)];
        let loads = SpatialMotionLoads {
            held: held(&[0.], &[ZERO; 2]),
            world_points: &[],
            connections: &springs,
        };
        let mut actual = state(vec![0.1], vec![-0.4]);
        let mut oracle = (0.1, -0.4);
        let initial_energy = 0.5 * mu * oracle.1 * oracle.1
            + 0.5 * springs[0].stiffness_n_per_m * oracle.0 * oracle.0;
        let mut loss = 0.;
        let mut max_evaluations = 0;
        for _ in 0..60 {
            let (q, v, mid) = midpoint_mode(
                oracle.0,
                oracle.1,
                springs[0].stiffness_n_per_m / mu,
                damping / mu,
                h,
            );
            let result = step_free_floating_spatial_implicit(
                &model,
                &actual,
                loads,
                h,
                ImplicitSpatialConfig {
                    residual_tolerance: 1e-8,
                    ..Default::default()
                },
                || false,
            )
            .unwrap();
            max_evaluations = max_evaluations.max(result.residual_evaluations);
            let s = &result.spatial;
            close(s.motion.state_after.joint_position[0], q, 2e-11);
            close(s.motion.state_after.joint_velocity[0], v, 4e-9);
            close(s.dissipated_energy, h * damping * mid * mid, 2e-7);
            close(
                s.motion.state_after.base.twist_body.linear.x,
                (-0.8 - 2. * v) / 5.,
                2e-9,
            );
            balance(&s.motion, 2e-11);
            assert!(result.scaled_residual_norm <= 1e-8);
            assert!(s.mechanical_work_residual.abs() < 2e-7);
            loss += s.dissipated_energy;
            actual = s.motion.state_after.clone();
            oracle = (q, v);
        }
        close(
            0.5 * mu * oracle.1 * oracle.1
                + 0.5 * springs[0].stiffness_n_per_m * oracle.0 * oracle.0
                + loss,
            initial_energy,
            3e-6,
        );
        eprintln!(
            "stiff mode h*omega=10, damping={damping}, max residual evaluations={max_evaluations}, loss={loss:e}"
        );
        // A stable midpoint response retains unresolved oscillation; no L-stable claim.
        if damping == 0. {
            assert!((oracle.0 * oracle.0 + (oracle.1 / 1000.).powi(2)).sqrt() > 0.099);
        }
    }
}

#[test]
fn resolved_implicit_oscillator_converges_second_order_to_independent_continuous_solution() {
    let model = slider();
    let springs = [connection(0, 1, 120., 6.)];
    let loads = SpatialMotionLoads {
        held: held(&[0.], &[ZERO; 2]),
        world_points: &[],
        connections: &springs,
    };
    let mu = 1.2;
    let alpha = springs[0].damping_ns_per_m / (2. * mu);
    let omega = (springs[0].stiffness_n_per_m / mu - alpha * alpha).sqrt();
    let t = 0.6;
    let a = 0.25;
    let b = (-0.15 + alpha * a) / omega;
    let exact_q = (-alpha * t).exp() * (a * (omega * t).cos() + b * (omega * t).sin());
    let exact_v = -alpha * exact_q
        + (-alpha * t).exp() * omega * (-a * (omega * t).sin() + b * (omega * t).cos());
    let mut errors = Vec::new();
    for steps in [30, 60, 120] {
        let mut current = state(vec![a], vec![-0.15]);
        for _ in 0..steps {
            let result = step_free_floating_spatial_implicit(
                &model,
                &current,
                loads,
                t / f64::from(steps),
                ImplicitSpatialConfig::default(),
                || false,
            )
            .unwrap();
            balance(&result.spatial.motion, 1e-12);
            assert!(result.spatial.mechanical_work_residual.abs() < 1e-10);
            current = result.spatial.motion.state_after;
        }
        errors.push(
            ((current.joint_position[0] - exact_q).powi(2)
                + (current.joint_velocity[0] - exact_v).powi(2))
            .sqrt(),
        );
    }
    eprintln!("implicit exact oscillator state errors: {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.9);
    }
}

#[test]
fn shared_link_connections_solve_both_independent_three_mass_modes_simultaneously() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "left",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(1., Vec3::default()),
        ),
        Link::new(
            "middle",
            Some(0),
            pose(Vec3::new(1., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(1., Vec3::default()),
        ),
        Link::new(
            "right",
            Some(0),
            pose(Vec3::new(2., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(1., Vec3::default()),
        ),
    ])
    .unwrap();
    let springs = [connection(0, 1, 2e4, 20.), connection(1, 2, 2e4, 20.)];
    let loads = SpatialMotionLoads {
        held: held(&[0.; 2], &[ZERO; 3]),
        world_points: &[],
        connections: &springs,
    };
    let a = 0.04;
    let b = -0.02;
    let va = -0.1;
    let vb = 0.08;
    let bulk = 2.3;
    // Physical mode shapes [1,0,-1] and [1,-2,1] have eigenvalues 1 and 3.
    let displacement = [a + b, -2. * b, -a + b];
    let velocity = [va + vb, -2. * vb, -va + vb];
    let mut original = state(
        vec![
            displacement[1] - displacement[0],
            displacement[2] - displacement[0],
        ],
        vec![velocity[1] - velocity[0], velocity[2] - velocity[0]],
    );
    original.base.world_from_base = pose(Vec3::new(displacement[0], 0., 0.));
    original.base.twist_body.linear = Vec3::new(bulk + velocity[0], 0., 0.);
    let h = 0.02;
    let result = step_free_floating_spatial_implicit(
        &model,
        &original,
        loads,
        h,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let (an, van, _) = midpoint_mode(a, va, 2e4, 20., h);
    let (bn, vbn, _) = midpoint_mode(b, vb, 6e4, 60., h);
    let q = [an + bn, -2. * bn, -an + bn];
    let v = [van + vbn, -2. * vbn, -van + vbn];
    let next = &result.spatial.motion.state_after;
    close(
        next.base.world_from_base.translation().x,
        q[0] + bulk * h,
        4e-13,
    );
    close(next.base.twist_body.linear.x, bulk + v[0], 2e-11);
    for j in 0..2 {
        close(next.joint_position[j], q[j + 1] - q[0], 4e-13);
        close(next.joint_velocity[j], v[j + 1] - v[0], 2e-11);
    }
    balance(&result.spatial.motion, 1e-12);
    assert!(result.spatial.mechanical_work_residual.abs() < 1e-10);
    let reversed = [springs[1], springs[0]];
    let again = step_free_floating_spatial_implicit(
        &model,
        &original,
        SpatialMotionLoads {
            connections: &reversed,
            ..loads
        },
        h,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    for (x, y) in next
        .joint_velocity
        .iter()
        .zip(&again.spatial.motion.state_after.joint_velocity)
    {
        close(*x, *y, 2e-11);
    }
}

fn mixed() -> (
    ArticulatedModel,
    FreeFloatingMotionState,
    [AxialConnection; 2],
) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3., Vec3::new(0.1, -0.2, 0.1)),
        ),
        Link::new(
            "hinge",
            Some(0),
            Se3::from_parts(rot(Vec3::new(0.2, 0.1, -0.1)), Vec3::new(0.8, 0.3, -0.2)).unwrap(),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            inertia(2., Vec3::new(0.2, 0.1, -0.1)),
        ),
        Link::new(
            "slide",
            Some(1),
            Se3::from_parts(rot(Vec3::new(-0.1, 0.2, 0.1)), Vec3::new(0.6, -0.2, 0.4)).unwrap(),
            JointModel::prismatic(Vec3::new(1., 0.3, 0.2), None).unwrap(),
            inertia(1., Vec3::new(0.1, -0.1, 0.2)),
        ),
    ])
    .unwrap();
    let mut current = state(vec![0.3, 0.2], vec![0.7, -0.15]);
    current.base = FreeFloatingBaseState::new(
        Se3::from_parts(rot(Vec3::new(0.3, -0.2, 0.1)), Vec3::new(-0.2, 0.3, 0.4)).unwrap(),
        Twist::new(Vec3::new(0.2, -0.3, 0.1), Vec3::new(0.4, -0.2, 0.1)),
    );
    let c = AxialConnection {
        a: point(0, Vec3::new(-0.1, 0.3, -0.2)),
        b: point(2, Vec3::new(0.2, -0.3, 0.1)),
        rest_length_m: 1.2,
        stiffness_n_per_m: 70.,
        damping_ns_per_m: 2.,
    };
    let d = AxialConnection {
        a: point(1, Vec3::new(0.3, 0.1, 0.2)),
        b: point(2, Vec3::new(0.2, -0.3, 0.1)),
        rest_length_m: 0.7,
        stiffness_n_per_m: 40.,
        damping_ns_per_m: 1.,
    };
    (model, current, [c, d])
}
fn state_error(a: &FreeFloatingMotionState, b: &FreeFloatingMotionState) -> f64 {
    let mut squared = (a.base.world_from_base.translation() - b.base.world_from_base.translation())
        .norm()
        .powi(2)
        + (a.base.twist_body.angular - b.base.twist_body.angular)
            .norm()
            .powi(2)
        + (a.base.twist_body.linear - b.base.twist_body.linear)
            .norm()
            .powi(2);
    for (x, y) in a
        .joint_position
        .iter()
        .chain(&a.joint_velocity)
        .zip(b.joint_position.iter().chain(&b.joint_velocity))
    {
        squared += (x - y) * (x - y);
    }
    for (x, y) in a
        .base
        .world_from_base
        .rotation()
        .matrix()
        .m
        .iter()
        .zip(b.base.world_from_base.rotation().matrix().m)
    {
        squared += (x - y) * (x - y);
    }
    squared.sqrt()
}

#[test]
fn changing_geometry_and_spatial_loads_converge_with_full_tree_momentum_and_work() {
    let (model, start, springs) = mixed();
    let applied = [WorldPointLoad {
        point: point(2, Vec3::new(0.2, -0.1, 0.3)),
        force_world: Vec3::new(1.2, -0.8, 0.6),
        couple_world: Vec3::new(0.1, -0.2, 0.3),
    }];
    let body_wrench = [
        ZERO,
        Wrench::new(Vec3::new(-0.15, 0.1, 0.05), Vec3::new(0.3, -0.2, 0.1)),
        ZERO,
    ];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: Vec3::new(0., -0.3, 0.1),
            ..held(&[0.2, -0.1], &body_wrench)
        },
        world_points: &applied,
        connections: &springs,
    };
    let t = 0.16;
    // A separately formulated explicit RK midpoint trajectory on a fine rung.
    let mut reference = start.clone();
    for _ in 0..4096 {
        reference = step_free_floating_spatial(&model, &reference, loads, t / 4096., || false)
            .unwrap()
            .motion
            .state_after;
    }
    let mut errors = Vec::new();
    let mut defects = Vec::new();
    for steps in [8, 16, 32] {
        let mut current = start.clone();
        let mut defect = 0.;
        for _ in 0..steps {
            let result = step_free_floating_spatial_implicit(
                &model,
                &current,
                loads,
                t / f64::from(steps),
                ImplicitSpatialConfig::default(),
                || false,
            )
            .unwrap();
            balance(&result.spatial.motion, 3e-12);
            assert!(result.spatial.dissipated_energy >= 0.);
            defect += result.spatial.mechanical_work_residual;
            current = result.spatial.motion.state_after;
        }
        errors.push(state_error(&current, &reference));
        defects.push(defect.abs());
    }
    eprintln!(
        "changing geometry implicit errors: {errors:?}; accumulated work defects: {defects:?}"
    );
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8);
    }
    for pair in defects.windows(2) {
        assert!(pair[0] / pair[1] > 3.5);
    }
}

#[test]
fn frame_changes_rate_scaling_and_large_common_motion_preserve_the_attachment_response() {
    let (model, original, springs) = mixed();
    let loads = SpatialMotionLoads {
        held: held(&[0.; 2], &[ZERO; 3]),
        world_points: &[],
        connections: &springs,
    };
    let a = step_free_floating_spatial_implicit(
        &model,
        &original,
        loads,
        0.01,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let transform =
        Se3::from_parts(rot(Vec3::new(-0.5, 0.6, 0.2)), Vec3::new(2., -3., 1.)).unwrap();
    let mut transformed = original.clone();
    transformed.base.world_from_base = transform.compose(original.base.world_from_base).unwrap();
    let b = step_free_floating_spatial_implicit(
        &model,
        &transformed,
        loads,
        0.01,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let mut expected = a.spatial.motion.state_after.clone();
    expected.base.world_from_base = transform.compose(expected.base.world_from_base).unwrap();
    assert!(state_error(&b.spatial.motion.state_after, &expected) < 2e-10);
    close(
        a.spatial.dissipated_energy,
        b.spatial.dissipated_energy,
        2e-11,
    );
    let scaled = step_free_floating_spatial_implicit(
        &model,
        &original,
        loads,
        0.01,
        ImplicitSpatialConfig {
            linear_velocity_scale: 0.1,
            angular_velocity_scale: 10.,
            residual_tolerance: 1e-11,
            ..Default::default()
        },
        || false,
    )
    .unwrap();
    assert!(
        state_error(
            &a.spatial.motion.state_after,
            &scaled.spatial.motion.state_after
        ) < 2e-10
    );
    let model = slider();
    let spring = [connection(0, 1, 1200., 6.)];
    let loads = SpatialMotionLoads {
        held: held(&[0.], &[ZERO; 2]),
        world_points: &[],
        connections: &spring,
    };
    let original = state(vec![0.25], vec![-0.15]);
    let mut fast = original.clone();
    fast.base.twist_body.linear.x = 1e9;
    let a = step_free_floating_spatial_implicit(
        &model,
        &original,
        loads,
        0.01,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let b = step_free_floating_spatial_implicit(
        &model,
        &fast,
        loads,
        0.01,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        a.spatial.motion.state_after.joint_velocity,
        b.spatial.motion.state_after.joint_velocity
    );
    assert_eq!(a.spatial.dissipated_energy, b.spatial.dissipated_energy);
    assert_eq!(
        a.spatial.motion.internal_body_work,
        b.spatial.motion.internal_body_work
    );
}

#[test]
fn a_single_body_under_world_force_has_exact_linear_impulse_and_work() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2., Vec3::default()),
    )])
    .unwrap();
    let mut original = state(vec![], vec![]);
    original.base.twist_body.linear = Vec3::new(1., -2., 3.);
    let force = Vec3::new(3., -1., 2.);
    let loads = SpatialMotionLoads {
        held: held(&[], &[ZERO]),
        world_points: &[WorldPointLoad {
            point: point(0, Vec3::default()),
            force_world: force,
            couple_world: Vec3::default(),
        }],
        connections: &[],
    };
    let result = step_free_floating_spatial_implicit(
        &model,
        &original,
        loads,
        0.2,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    vector(
        result.spatial.motion.state_after.base.twist_body.linear,
        original.base.twist_body.linear + force.scale(0.1),
        1e-14,
    );
    vector(
        result
            .spatial
            .motion
            .state_after
            .base
            .world_from_base
            .translation(),
        original.base.twist_body.linear.scale(0.2) + force.scale(0.01),
        1e-14,
    );
    close(result.spatial.mechanical_work_residual, 0., 1e-13);
    assert_eq!(result.nonlinear.iterations, 0);
    assert_eq!(result.residual_evaluations, 2);
    balance(&result.spatial.motion, 2e-13);
}

#[test]
fn cancellation_budgets_bad_controls_and_late_geometry_refuse_atomically() {
    let (model, original, springs) = mixed();
    let before = original.clone();
    let loads = SpatialMotionLoads {
        held: held(&[0.; 2], &[ZERO; 3]),
        world_points: &[],
        connections: &springs,
    };
    let config = ImplicitSpatialConfig::default();
    let mut polls = 0;
    let expected =
        step_free_floating_spatial_implicit(&model, &original, loads, 0.04, config, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            step_free_floating_spatial_implicit(&model, &original, loads, 0.04, config, || {
                seen += 1;
                seen == stop
            }),
            Err(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::Cancelled
            ))
        ));
        assert_eq!(original, before);
    }
    assert_eq!(
        expected,
        step_free_floating_spatial_implicit(&model, &original, loads, 0.04, config, || false)
            .unwrap()
    );
    assert!(matches!(
        step_free_floating_spatial_implicit(
            &model,
            &original,
            loads,
            0.04,
            ImplicitSpatialConfig {
                maximum_evaluations: 2,
                ..config
            },
            || false
        ),
        Err(ImplicitSpatialError::EvaluationBudget {
            used: 2,
            maximum: 2
        })
    ));
    let late_budget = expected.residual_evaluations - 1;
    assert!(
        matches!(step_free_floating_spatial_implicit(&model,&original,loads,0.04,
        ImplicitSpatialConfig { maximum_evaluations:late_budget,..config },||false),
        Err(ImplicitSpatialError::EvaluationBudget { used,maximum }) if used==late_budget && maximum==late_budget)
    );
    assert!(matches!(
        step_free_floating_spatial_implicit(
            &model,
            &original,
            loads,
            0.04,
            ImplicitSpatialConfig {
                maximum_iterations: 1,
                ..config
            },
            || false
        ),
        Err(ImplicitSpatialError::Unresolved(_))
    ));
    for bad in [
        ImplicitSpatialConfig {
            angular_velocity_scale: f64::NAN,
            ..config
        },
        ImplicitSpatialConfig {
            linear_velocity_scale: 0.,
            ..config
        },
        ImplicitSpatialConfig {
            residual_tolerance: 0.,
            ..config
        },
        ImplicitSpatialConfig {
            difference_step: f64::INFINITY,
            ..config
        },
        ImplicitSpatialConfig {
            maximum_iterations: 65,
            ..config
        },
        ImplicitSpatialConfig {
            krylov_restart: 65,
            ..config
        },
        ImplicitSpatialConfig {
            maximum_evaluations: MAX_IMPLICIT_SPATIAL_EVALUATIONS + 1,
            ..config
        },
    ] {
        assert!(matches!(
            step_free_floating_spatial_implicit(&model, &original, loads, 0.04, bad, || false),
            Err(ImplicitSpatialError::Invalid(_))
        ));
    }
    let model = slider();
    let original = state(vec![0.], vec![-2.]);
    let spring = [connection(0, 1, 0., 0.)];
    let loads = SpatialMotionLoads {
        held: held(&[0.], &[ZERO; 2]),
        world_points: &[],
        connections: &spring,
    };
    assert!(matches!(
        step_free_floating_spatial_implicit(&model, &original, loads, 0.5, config, || false),
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::InvalidLoad(_)
        ))
    ));
    assert_eq!(original.joint_position, [0.]);
}
