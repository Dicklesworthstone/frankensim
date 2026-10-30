//! Physical rigid loop motion, independent reductions, and actual reaction work.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    DistanceLoop, DistanceLoopConfig, DistanceLoopError, DistanceLoopStep, MAX_DISTANCE_LOOPS,
    step_free_floating_distance_loops, step_free_floating_distance_loops_with_bushings,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialError, step_free_floating_spatial_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::{
    MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics, free_floating_momentum,
};

const ZERO: Wrench = Wrench::new(
    Vec3 {
        x: 0.,
        y: 0.,
        z: 0.,
    },
    Vec3 {
        x: 0.,
        y: 0.,
        z: 0.,
    },
);
fn close(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t,
        "{a:.16e} != {b:.16e}, error {}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, t: f64) {
    assert!(
        (a - b).norm() <= t,
        "{a:?} != {b:?}, error {}",
        (a - b).norm()
    );
}
fn rotation(w: Vec3) -> So3 {
    So3::exp(So3Tangent::new(w)).unwrap()
}
fn pose(v: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), v).unwrap()
}
fn inertia(m: f64) -> SpatialInertia {
    SpatialInertia::new(
        m,
        Vec3::default(),
        Mat3 {
            m: [2., 0.2, 0.1, 0.2, 2.5, 0.15, 0.1, 0.15, 3.],
        },
    )
    .unwrap()
}
fn point(link: usize) -> MaterialPoint {
    MaterialPoint {
        link,
        point_body: Vec3::default(),
    }
}
fn dist(a: usize, b: usize, l: f64) -> DistanceLoop {
    DistanceLoop {
        a: point(a),
        b: point(b),
        length_m: l,
    }
}
fn held<'a>(q: &'a [f64], w: &'a [Wrench]) -> FreeFloatingMotionLoads<'a> {
    FreeFloatingMotionLoads {
        generalized_force: q,
        external_body_wrench: w,
        gravity_world: Vec3::default(),
    }
}
fn still(q: Vec<f64>, v: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: q,
        joint_velocity: v,
    }
}
fn account(s: &DistanceLoopStep) {
    assert!(s.position_residual_after_m.iter().all(|x| x.abs() < 1e-9));
    assert!(s.velocity_residual_after_mps.iter().all(|x| x.abs() < 1e-9));
    assert!(s.linear_momentum_residual_world.norm() < 2e-11);
    assert!(s.angular_momentum_residual_world.norm() < 3e-11);
    close(
        s.mechanical_work_residual,
        s.constraint_impulse_work
            + s.midpoint.mechanical_work_residual
            + s.initial_reaction.kinetic_work_residual
            + s.final_reaction.kinetic_work_residual,
        2e-11,
    );
}
fn slider() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(Vec3::new(1., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2.),
        ),
    ])
    .unwrap()
}

#[test]
fn rigid_locked_slider_shares_an_external_force_with_exact_reaction_impulses() {
    let model = slider();
    let loops = [dist(0, 1, 1.)];
    let mut state = still(vec![0.], vec![0.]);
    state.base.twist_body.linear = Vec3::new(0.7, 0., 0.);
    let force = 3.2;
    let h = 0.04;
    let load = [WorldPointLoad {
        point: point(1),
        force_world: Vec3::new(force, 0., 0.),
        couple_world: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: held(&[0.], &[ZERO; 2]),
        world_points: &load,
        connections: &[],
    };
    let result = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &loops,
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&result);
    close(result.state_after.joint_position[0], 0., 2e-13);
    close(result.state_after.joint_velocity[0], 0., 2e-13);
    close(
        result.state_after.base.twist_body.linear.x,
        0.7 + force * h / 5.,
        2e-13,
    );
    close(
        result.state_after.base.world_from_base.translation().x,
        0.7 * h + 0.5 * force * h * h / 5.,
        2e-13,
    );
    close(result.initial_impulse_ns[0], -0.3 * force * h, 2e-12);
    close(result.final_impulse_ns[0], -0.3 * force * h, 2e-12);
    assert!(result.initial_reaction.impulse_work > 0.);
    assert!(result.final_reaction.impulse_work < 0.);
    close(result.constraint_impulse_work, 0., 2e-13);
    close(result.mechanical_work_residual, 0., 2e-13);
}

fn triangle() -> (
    ArticulatedModel,
    ArticulatedModel,
    FreeFloatingMotionState,
    FreeFloatingMotionState,
    DistanceLoop,
) {
    let root = Link::new(
        "base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(3.),
    );
    let rod = SpatialInertia::new(
        1.,
        Vec3::new(0.5, 0., 0.),
        Mat3 {
            m: [0.1, 0., 0., 0., 0.2, 0., 0., 0., 0.25],
        },
    )
    .unwrap();
    let hinge = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
    let model = ArticulatedModel::new(vec![
        root.clone(),
        Link::new("arm_a", Some(0), Se3::identity(), hinge, rod),
        Link::new("arm_b", Some(0), Se3::identity(), hinge, rod),
    ])
    .unwrap();
    // Independent parallel-axis sum of two original arms at relative pi/2.
    let composite = SpatialInertia::new(
        2.,
        Vec3::new(0.25, 0.25, 0.),
        Mat3 {
            m: [0.425, 0.125, 0., 0.125, 0.425, 0., 0., 0., 0.75],
        },
    )
    .unwrap();
    let reduced = ArticulatedModel::new(vec![
        root,
        Link::new(
            "complete_triangle",
            Some(0),
            Se3::identity(),
            hinge,
            composite,
        ),
    ])
    .unwrap();
    let q = 0.4;
    let v = 0.7;
    let half_pi = core::f64::consts::FRAC_PI_2;
    let mut state = still(vec![q, q + half_pi], vec![v, v]);
    state.base = FreeFloatingBaseState::new(
        Se3::from_parts(
            rotation(Vec3::new(0.3, -0.2, 0.1)),
            Vec3::new(0.2, -0.3, 0.4),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.2, 0.1, -0.1)),
    );
    let reference = FreeFloatingMotionState {
        joint_position: vec![q],
        joint_velocity: vec![v],
        base: state.base,
    };
    let constraint = DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::new(1., 0., 0.),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::new(1., 0., 0.),
        },
        length_m: 2.0_f64.sqrt(),
    };
    (model, reduced, state, reference, constraint)
}
fn reduced_error(a: &FreeFloatingMotionState, b: &FreeFloatingMotionState) -> f64 {
    let mut error = (a.base.world_from_base.translation() - b.base.world_from_base.translation())
        .norm()
        .powi(2)
        + (a.base.twist_body.linear - b.base.twist_body.linear)
            .norm()
            .powi(2)
        + (a.base.twist_body.angular - b.base.twist_body.angular)
            .norm()
            .powi(2)
        + (a.joint_position[0] - b.joint_position[0]).powi(2)
        + (a.joint_velocity[0] - b.joint_velocity[0]).powi(2);
    for (x, y) in a
        .base
        .world_from_base
        .rotation()
        .matrix()
        .m
        .iter()
        .zip(b.base.world_from_base.rotation().matrix().m)
    {
        error += (x - y).powi(2);
    }
    error.sqrt()
}

#[test]
fn a_loaded_closed_triangle_converges_to_its_independently_reduced_full_inertia_body() {
    let (model, reduced, start, mut reference, constraint) = triangle();
    let loops = [constraint];
    let g = Vec3::new(0., -0.4, 0.);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: g,
            ..held(&[0.8, -0.2], &[ZERO; 3])
        },
        world_points: &[],
        connections: &[],
    };
    let reduced_loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: g,
            ..held(&[0.6], &[ZERO; 2])
        },
        world_points: &[],
        connections: &[],
    };
    let original = free_floating_momentum(
        &model,
        start.base,
        &start.joint_position,
        &start.joint_velocity,
    )
    .unwrap();
    let expected = free_floating_momentum(
        &reduced,
        reference.base,
        &reference.joint_position,
        &reference.joint_velocity,
    )
    .unwrap();
    vector(
        original.linear_momentum_world,
        expected.linear_momentum_world,
        2e-13,
    );
    vector(
        original.angular_momentum_world,
        expected.angular_momentum_world,
        2e-13,
    );
    close(original.kinetic_energy, expected.kinetic_energy, 2e-13);
    let time = 0.4;
    for _ in 0..1024 {
        reference = step_free_floating_spatial_implicit(
            &reduced,
            &reference,
            reduced_loads,
            time / 1024.,
            Default::default(),
            || false,
        )
        .unwrap()
        .spatial
        .motion
        .state_after;
    }
    let mut errors = Vec::new();
    let mut defects = Vec::new();
    for count in [12, 24, 48] {
        let mut state = start.clone();
        let mut defect = 0.;
        for _ in 0..count {
            let step = step_free_floating_distance_loops(
                &model,
                &state,
                loads,
                &loops,
                time / f64::from(count),
                Default::default(),
                || false,
            )
            .unwrap();
            account(&step);
            close(
                step.state_after.joint_position[1] - step.state_after.joint_position[0],
                core::f64::consts::FRAC_PI_2,
                2e-10,
            );
            close(
                step.state_after.joint_velocity[1],
                step.state_after.joint_velocity[0],
                2e-10,
            );
            defect += step.mechanical_work_residual;
            state = step.state_after;
        }
        errors.push(reduced_error(&state, &reference));
        defects.push(defect.abs());
    }
    eprintln!("rigid triangle reduced-body errors {errors:?}; mechanical work defects {defects:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.5, "{errors:?}");
    }
}

fn circle() -> (ArticulatedModel, FreeFloatingMotionState) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.),
        ),
        Link::new(
            "x_slide",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(1.),
        ),
        Link::new(
            "y_slide",
            Some(1),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0., 1., 0.), None).unwrap(),
            inertia(1.),
        ),
    ])
    .unwrap();
    let mut state = still(vec![0.8, 0.6], vec![-0.6, 0.8]);
    state.base.twist_body = Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.4, -0.3, 0.2));
    (model, state)
}

#[test]
fn nonlinear_loop_keeps_actual_world_material_geometry_and_tangent_velocities() {
    let (model, mut state) = circle();
    let loops = [dist(0, 2, 1.)];
    let point_load = [WorldPointLoad {
        point: point(2),
        force_world: Vec3::new(0.5, -0.2, 0.1),
        couple_world: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: Vec3::new(0., -0.4, 0.),
            ..held(&[0.2, -0.1], &[ZERO; 3])
        },
        world_points: &point_load,
        connections: &[],
    };
    let mut max_gap = 0.0_f64;
    let mut max_rate = 0.0_f64;
    let mut reaction_work = 0.;
    for _ in 0..40 {
        let step = step_free_floating_distance_loops(
            &model,
            &state,
            loads,
            &loops,
            0.01,
            Default::default(),
            || false,
        )
        .unwrap();
        account(&step);
        assert!(step.initial_reaction.impulse_work >= -1e-11);
        assert!(step.final_reaction.impulse_work <= 1e-11);
        state = step.state_after;
        let shape = forward_kinematics(
            &model,
            BaseState::prescribed(
                state.base.world_from_base,
                state.base.twist_body,
                Twist::zero(),
            ),
            &state.joint_position,
            &state.joint_velocity,
        )
        .unwrap();
        let delta = shape.world_from_link[2].translation() - shape.world_from_link[0].translation();
        let va = shape.world_from_link[0]
            .rotation()
            .rotate(shape.body_twist[0].linear)
            .unwrap();
        let vb = shape.world_from_link[2]
            .rotation()
            .rotate(shape.body_twist[2].linear)
            .unwrap();
        let gap = delta.norm() - 1.;
        let rate = delta.dot(vb - va) / delta.norm();
        max_gap = max_gap.max(gap.abs());
        max_rate = max_rate.max(rate.abs());
        close(gap, step.position_residual_after_m[0], 2e-13);
        close(rate, step.velocity_residual_after_mps[0], 2e-13);
        reaction_work += step.constraint_impulse_work;
    }
    assert!(state.joint_position[0] < 0.7);
    assert!(max_gap < 1e-10 && max_rate < 1e-10);
    eprintln!(
        "nonlinear circle maximum gap={max_gap:e}, rate={max_rate:e}, finite constraint work={reaction_work:e}"
    );
}

#[test]
fn simultaneous_shared_link_loops_keep_the_exact_force_sharing_and_permutation() {
    let model = ArticulatedModel::new(vec![
        Link::new("a", None, Se3::identity(), JointModel::FIXED, inertia(3.)),
        Link::new(
            "b",
            Some(0),
            pose(Vec3::new(1., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2.),
        ),
        Link::new(
            "c",
            Some(0),
            pose(Vec3::new(2., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(1.),
        ),
    ])
    .unwrap();
    let loops = [dist(0, 1, 1.), dist(1, 2, 1.)];
    let state = still(vec![0.; 2], vec![0.; 2]);
    let f = 7.;
    let h = 0.02;
    let applied = [WorldPointLoad {
        point: point(2),
        force_world: Vec3::new(f, 0., 0.),
        couple_world: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: held(&[0.; 2], &[ZERO; 3]),
        world_points: &applied,
        connections: &[],
    };
    let a = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &loops,
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    let b = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &[loops[1], loops[0]],
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&a);
    account(&b);
    close(a.state_after.base.twist_body.linear.x, f * h / 6., 3e-13);
    for q in &a.state_after.joint_position {
        close(*q, 0., 2e-13);
    }
    for v in &a.state_after.joint_velocity {
        close(*v, 0., 2e-13);
    }
    close(
        a.initial_impulse_ns[0] + a.final_impulse_ns[0],
        -0.5 * f * h,
        2e-12,
    );
    close(
        a.initial_impulse_ns[1] + a.final_impulse_ns[1],
        -(5. / 6.) * f * h,
        2e-12,
    );
    vector(
        a.state_after.base.twist_body.linear,
        b.state_after.base.twist_body.linear,
        2e-12,
    );
    close(a.initial_impulse_ns[0], b.initial_impulse_ns[1], 2e-12);
    close(a.initial_impulse_ns[1], b.initial_impulse_ns[0], 2e-12);
}

#[test]
fn rigid_support_and_finite_bushing_storage_damping_share_one_original_tree() {
    let (model, mut state) = circle();
    let loops = [dist(0, 2, 1.)];
    let bushings = [SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::new(1., 0., 0.),
        translation_stiffness: Vec3::new(20., 40., 0.),
        translation_damping: Vec3::new(0.5, 0.8, 0.),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: held(&[0.; 2], &[ZERO; 3]),
        world_points: &[],
        connections: &[],
    };
    let plain = step_free_floating_distance_loops(
        &model,
        &state,
        loads,
        &loops,
        0.005,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        plain,
        step_free_floating_distance_loops_with_bushings(
            &model,
            &state,
            loads,
            &loops,
            &[],
            0.005,
            Default::default(),
            || false
        )
        .unwrap()
    );
    let initial = free_floating_momentum(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let initial_energy = initial.kinetic_energy + 7.6;
    let mut loss = 0.;
    let mut defect = 0.;
    for step in 0..24 {
        let result = step_free_floating_distance_loops_with_bushings(
            &model,
            &state,
            loads,
            &loops,
            &bushings,
            0.005,
            Default::default(),
            || false,
        )
        .unwrap();
        account(&result);
        if step == 0 {
            close(result.midpoint.stored_energy_before, 7.6, 1e-13);
            assert!(
                (result.state_after.joint_velocity[0] - plain.state_after.joint_velocity[0]).abs()
                    > 0.005
            );
        }
        loss += result.midpoint.dissipated_energy;
        defect += result.mechanical_work_residual;
        assert!(result.midpoint.dissipated_energy > 0.);
        close(
            result.final_reaction.momentum_after.kinetic_energy
                + result.midpoint.stored_energy_after
                + loss,
            initial_energy + defect,
            2e-11,
        );
        state = result.state_after;
    }
    assert!(loss > 0.01);
    eprintln!(
        "rigid circle with anisotropic finite bushing: damping loss={loss:e}, mechanical defect={defect:e}"
    );
}

#[test]
fn rigid_loop_reaction_and_motion_are_scene_covariant() {
    let (model, start) = circle();
    let loops = [dist(0, 2, 1.)];
    let g = Vec3::new(0., -0.3, 0.);
    let force = Vec3::new(0.5, -0.2, 0.1);
    let applied = [WorldPointLoad {
        point: point(2),
        force_world: force,
        couple_world: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: g,
            ..held(&[0.2, -0.1], &[ZERO; 3])
        },
        world_points: &applied,
        connections: &[],
    };
    let a = step_free_floating_distance_loops(
        &model,
        &start,
        loads,
        &loops,
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    let transform =
        Se3::from_parts(rotation(Vec3::new(-0.3, 0.2, 0.4)), Vec3::new(2., -1., 0.3)).unwrap();
    let mut transformed = start.clone();
    transformed.base.world_from_base = transform.compose(start.base.world_from_base).unwrap();
    let rotated = [WorldPointLoad {
        force_world: transform.rotation().rotate(force).unwrap(),
        ..applied[0]
    }];
    let b = step_free_floating_distance_loops(
        &model,
        &transformed,
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                gravity_world: transform.rotation().rotate(g).unwrap(),
                ..loads.held
            },
            world_points: &rotated,
            ..loads
        },
        &loops,
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    let expected = transform
        .compose(a.state_after.base.world_from_base)
        .unwrap();
    vector(
        b.state_after.base.world_from_base.translation(),
        expected.translation(),
        2e-11,
    );
    vector(
        a.state_after.base.twist_body.angular,
        b.state_after.base.twist_body.angular,
        2e-11,
    );
    vector(
        a.state_after.base.twist_body.linear,
        b.state_after.base.twist_body.linear,
        2e-11,
    );
    for (x, y) in a
        .state_after
        .joint_velocity
        .iter()
        .zip(&b.state_after.joint_velocity)
    {
        close(*x, *y, 2e-11);
    }
    close(a.initial_impulse_ns[0], b.initial_impulse_ns[0], 2e-11);
    close(a.constraint_impulse_work, b.constraint_impulse_work, 2e-12);
}

#[test]
fn inconsistent_redundant_and_already_imposed_distance_constraints_refuse() {
    let model = slider();
    let state = still(vec![0.], vec![0.]);
    let constraint = dist(0, 1, 1.);
    let loads = SpatialMotionLoads {
        held: held(&[0.], &[ZERO; 2]),
        world_points: &[],
        connections: &[],
    };
    let run = |state: &FreeFloatingMotionState, c: &[DistanceLoop]| {
        step_free_floating_distance_loops(&model, state, loads, c, 0.01, Default::default(), || {
            false
        })
    };
    assert!(matches!(
        run(
            &state,
            &[DistanceLoop {
                length_m: 1.01,
                ..constraint
            }]
        ),
        Err(DistanceLoopError::Position {
            stage: "initial",
            ..
        })
    ));
    assert!(matches!(
        run(&still(vec![0.], vec![0.1]), &[constraint]),
        Err(DistanceLoopError::Velocity {
            stage: "initial",
            ..
        })
    ));
    assert!(matches!(
        run(&state, &[constraint, constraint]),
        Err(DistanceLoopError::Mobility {
            stage: "initial",
            ..
        })
    ));
    assert!(matches!(
        run(&state, &[]),
        Err(DistanceLoopError::Invalid(_))
    ));
    assert!(matches!(
        run(&state, &vec![constraint; MAX_DISTANCE_LOOPS + 1]),
        Err(DistanceLoopError::Invalid(_))
    ));
    let model = ArticulatedModel::new(vec![
        Link::new(
            "root",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.),
        ),
        Link::new(
            "hinged",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            inertia(1.),
        ),
    ])
    .unwrap();
    let redundant = DistanceLoop {
        a: point(0),
        b: MaterialPoint {
            link: 1,
            point_body: Vec3::new(1., 0., 0.),
        },
        length_m: 1.,
    };
    assert!(matches!(
        step_free_floating_distance_loops(
            &model,
            &state,
            loads,
            &[redundant],
            0.01,
            Default::default(),
            || false
        ),
        Err(DistanceLoopError::Mobility {
            stage: "initial",
            ..
        })
    ));
}

#[test]
fn late_cancellation_final_solve_budget_and_bad_geometry_keep_inputs_unchanged() {
    let (model, state) = circle();
    let original = state.clone();
    let loops = [dist(0, 2, 1.)];
    let loads = SpatialMotionLoads {
        held: held(&[0.2, -0.1], &[ZERO; 3]),
        world_points: &[],
        connections: &[],
    };
    let config = DistanceLoopConfig::default();
    let mut polls = 0;
    let expected =
        step_free_floating_distance_loops(&model, &state, loads, &loops, 0.02, config, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            step_free_floating_distance_loops(&model, &state, loads, &loops, 0.02, config, || {
                seen += 1;
                seen == stop
            }),
            Err(DistanceLoopError::Motion(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::Cancelled
            )))
        ));
        assert_eq!(state, original);
    }
    let budget = expected.residual_evaluations - 1;
    assert!(matches!(
        step_free_floating_distance_loops(
            &model,
            &state,
            loads,
            &loops,
            0.02,
            DistanceLoopConfig {
                nonlinear:
                    fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialConfig {
                        maximum_evaluations: budget,
                        ..config.nonlinear
                    },
                ..config
            },
            || false
        ),
        Err(DistanceLoopError::Motion(
            ImplicitSpatialError::EvaluationBudget { .. }
        ))
    ));
    assert_eq!(
        expected,
        step_free_floating_distance_loops(&model, &state, loads, &loops, 0.02, config, || false)
            .unwrap()
    );
    for bad in [
        DistanceLoopConfig {
            position_tolerance_m: f64::NAN,
            ..config
        },
        DistanceLoopConfig {
            velocity_tolerance_mps: 0.,
            ..config
        },
        DistanceLoopConfig {
            minimum_scaled_pivot: 1.,
            ..config
        },
        DistanceLoopConfig {
            minimum_scaled_pivot: f64::EPSILON,
            ..config
        },
    ] {
        assert!(matches!(
            step_free_floating_distance_loops(&model, &state, loads, &loops, 0.02, bad, || false),
            Err(DistanceLoopError::Invalid(_))
        ));
    }
    assert!(matches!(
        step_free_floating_distance_loops(
            &model,
            &state,
            loads,
            &[DistanceLoop {
                length_m: 0.,
                ..loops[0]
            }],
            0.02,
            config,
            || false
        ),
        Err(DistanceLoopError::Invalid(_))
    ));
}
