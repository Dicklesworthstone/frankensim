//! Physical bushing forces, finite-rotation work, native recoil and stiff modes.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MAX_BUSHING_ROTATION_RADIANS, MaterialFrame, SpatialBushing, evaluate_bushing_forces,
    step_free_floating_bushings, step_free_floating_bushings_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, step_free_floating_spatial_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad, step_free_floating_spatial,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics,
};

const Z: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};
const ZERO: Wrench = Wrench {
    torque: Z,
    force: Z,
};
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; error {:.4e}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}
fn rot(r: Vec3) -> So3 {
    So3::exp(So3Tangent::new(r)).unwrap()
}
fn pose(r: Vec3, x: Vec3) -> Se3 {
    Se3::from_parts(rot(r), x).unwrap()
}
fn inertia(m: f64, i: f64, com: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        m,
        com,
        Mat3 {
            m: [i, 0., 0., 0., i, 0., 0., 0., i],
        },
    )
    .unwrap()
}
fn frame(link: usize, r: Vec3, x: Vec3) -> MaterialFrame {
    MaterialFrame {
        link,
        link_from_attachment: pose(r, x),
    }
}
fn element(a: MaterialFrame, b: MaterialFrame) -> SpatialBushing {
    SpatialBushing {
        a,
        b,
        rest_translation_a: Z,
        translation_stiffness: Z,
        translation_damping: Z,
        rotation_stiffness: Z,
        rotation_damping: Z,
    }
}
fn held<'a>(force: &'a [f64], wrench: &'a [Wrench]) -> SpatialMotionLoads<'a> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: force,
            external_body_wrench: wrench,
            gravity_world: Z,
        },
        world_points: &[],
        connections: &[],
    }
}
fn balance(s: &FreeFloatingMotionStep, eps: f64) {
    vector(
        s.momentum_after.linear_momentum_world - s.momentum_before.linear_momentum_world,
        s.external_linear_impulse_world,
        eps,
    );
    vector(
        s.momentum_after.angular_momentum_world - s.momentum_before.angular_momentum_world,
        s.external_angular_impulse_world,
        eps,
    );
}
fn mixed() -> (
    ArticulatedModel,
    FreeFloatingMotionState,
    [SpatialBushing; 2],
) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3., 2., Vec3::new(0.1, -0.2, 0.1)),
        ),
        Link::new(
            "hinge",
            Some(0),
            pose(Vec3::new(0.2, 0.1, -0.1), Vec3::new(0.8, 0.3, -0.2)),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            inertia(2., 1.5, Vec3::new(0.2, 0.1, -0.1)),
        ),
        Link::new(
            "slide",
            Some(1),
            pose(Vec3::new(-0.1, 0.2, 0.1), Vec3::new(0.6, -0.2, 0.4)),
            JointModel::prismatic(Vec3::new(1., 0.3, 0.2), None).unwrap(),
            inertia(1., 1., Vec3::new(0.1, -0.1, 0.2)),
        ),
        Link::new(
            "branch",
            Some(0),
            pose(Vec3::new(0.1, -0.1, 0.2), Vec3::new(-0.3, 0.7, 0.2)),
            JointModel::revolute(Vec3::new(0.3, 1., 0.2), None).unwrap(),
            inertia(1.4, 1.2, Vec3::new(-0.2, 0.1, 0.2)),
        ),
    ])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(Vec3::new(0.3, -0.2, 0.1), Vec3::new(-0.2, 0.3, 0.4)),
            Twist::new(Vec3::new(0.2, -0.3, 0.1), Vec3::new(0.4, -0.2, 0.1)),
        ),
        joint_position: vec![0.7, 0.2, -0.4],
        joint_velocity: vec![0.7, -0.15, 0.3],
    };
    let a = SpatialBushing {
        rest_translation_a: Vec3::new(0.6, 0.1, -0.2),
        translation_stiffness: Vec3::new(20., 35., 50.),
        translation_damping: Vec3::new(0.6, 1.2, 0.3),
        rotation_stiffness: Vec3::new(12., 18., 7.),
        rotation_damping: Vec3::new(0.4, 0.7, 0.2),
        ..element(
            frame(0, Vec3::new(0.2, -0.1, 0.3), Vec3::new(-0.2, 0.1, 0.3)),
            frame(2, Vec3::new(-0.2, 0.4, -0.1), Vec3::new(0.3, -0.1, 0.2)),
        )
    };
    let b = SpatialBushing {
        rest_translation_a: Vec3::new(-0.4, 0.6, 0.2),
        translation_stiffness: Vec3::new(15., 10., 25.),
        translation_damping: Vec3::new(0.2, 0.7, 0.3),
        rotation_stiffness: Vec3::new(6., 3., 9.),
        rotation_damping: Vec3::new(0.2, 0.1, 0.4),
        ..element(
            frame(1, Vec3::new(0.1, 0.2, -0.2), Vec3::new(0.1, -0.2, 0.1)),
            frame(3, Vec3::new(-0.1, 0.3, 0.2), Vec3::new(0.2, 0.1, -0.1)),
        )
    };
    (model, state, [a, b])
}

fn reference_energy(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    bs: &[SpatialBushing],
) -> f64 {
    let fk = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    bs.iter()
        .map(|b| {
            let a = fk.world_from_link[b.a.link]
                .compose(b.a.link_from_attachment)
                .unwrap();
            let c = fk.world_from_link[b.b.link]
                .compose(b.b.link_from_attachment)
                .unwrap();
            let relative = a.inverse().unwrap().compose(c).unwrap();
            let q = relative.translation() - b.rest_translation_a;
            let r = relative.rotation().log().angular;
            let product = |v: Vec3, k: Vec3| k.x * v.x * v.x + k.y * v.y * v.y + k.z * v.z * v.z;
            0.5 * (product(q, b.translation_stiffness) + product(r, b.rotation_stiffness))
        })
        .sum()
}

#[test]
fn finite_anisotropic_forces_match_independent_energy_gradient_and_cancel_spatial_moment() {
    let (model, state, bs) = mixed();
    let response = evaluate_bushing_forces(&model, &state, &bs, || false).unwrap();
    close(
        response.stored_energy,
        reference_energy(&model, &state, &bs),
        2e-13,
    );
    let fk = forward_kinematics(
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
    let mut force = Z;
    let mut torque = Z;
    let mut power = 0.;
    for (i, w) in response.body_wrench.iter().enumerate() {
        let pose = fk.world_from_link[i];
        let f = pose.rotation().rotate(w.force).unwrap();
        force = force + f;
        torque = torque + pose.rotation().rotate(w.torque).unwrap() + pose.translation().cross(f);
        power += w.pairing(fk.body_twist[i]);
    }
    vector(force, Z, 2e-13);
    vector(torque, Z, 2e-13);
    close(power, response.internal_power, 2e-13);
    assert!(response.dissipation_power > 0.);
    // Independent finite difference of stored energy along the actual shape
    // velocity checks the geometric force couple AND finite-log torque map.
    let epsilon = 2e-6;
    let mut plus = state.clone();
    let mut minus = state.clone();
    for i in 0..model.dof_count() {
        plus.joint_position[i] += epsilon * state.joint_velocity[i];
        minus.joint_position[i] -= epsilon * state.joint_velocity[i];
    }
    let derivative = (reference_energy(&model, &plus, &bs) - reference_energy(&model, &minus, &bs))
        / (2. * epsilon);
    close(power + derivative + response.dissipation_power, 0., 3e-8);
    // Per-coordinate virtual work independently checks each generalized load.
    let conservative = bs.map(|b| SpatialBushing {
        translation_damping: Z,
        rotation_damping: Z,
        ..b
    });
    let response = evaluate_bushing_forces(&model, &state, &conservative, || false).unwrap();
    for j in 0..model.dof_count() {
        let mut velocity = vec![0.; model.dof_count()];
        velocity[j] = 1.;
        let fk = forward_kinematics(
            &model,
            BaseState::stationary(Se3::identity()),
            &state.joint_position,
            &velocity,
        )
        .unwrap();
        let load: f64 = response
            .body_wrench
            .iter()
            .zip(&fk.body_twist)
            .map(|(w, v)| w.pairing(*v))
            .sum();
        let mut plus = state.clone();
        let mut minus = state.clone();
        plus.joint_position[j] += epsilon;
        minus.joint_position[j] -= epsilon;
        let gradient = (reference_energy(&model, &plus, &bs)
            - reference_energy(&model, &minus, &bs))
            / (2. * epsilon);
        close(load, -gradient, 3e-8);
    }
}

#[test]
fn relative_law_is_exactly_independent_of_large_common_rigid_motion() {
    let (model, state, bs) = mixed();
    let first = evaluate_bushing_forces(&model, &state, &bs, || false).unwrap();
    let mut moving = state.clone();
    // Both poses satisfy native motor admission. Rotated dual quaternions at
    // enormous offsets can lose their Study condition before reaching this law.
    for world_pose in [
        pose(Vec3::new(-0.8, 0.5, 0.4), Vec3::new(1., -2., 3.)),
        pose(Z, Vec3::new(1e12, -2e12, 3e12)),
    ] {
        moving.base = FreeFloatingBaseState::new(
            world_pose,
            Twist::new(Vec3::new(1e8, 2e8, -3e8), Vec3::new(1e10, -2e10, 3e10)),
        );
        assert_eq!(
            first,
            evaluate_bushing_forces(&model, &moving, &bs, || false).unwrap()
        );
    }
}

#[test]
fn stiff_translation_and_coincident_torsion_match_two_body_reduced_mass_oscillators() {
    for torsion in [false, true] {
        for damping in [0., 240.] {
            let model = ArticulatedModel::new(vec![
                Link::new(
                    "base",
                    None,
                    Se3::identity(),
                    JointModel::FIXED,
                    inertia(3., 2., Z),
                ),
                Link::new(
                    "moving",
                    Some(0),
                    pose(Z, if torsion { Z } else { Vec3::new(1., 0., 0.) }),
                    if torsion {
                        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap()
                    } else {
                        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap()
                    },
                    inertia(2., 3., Z),
                ),
            ])
            .unwrap();
            let mu = 1.2;
            let k = 1.2e6;
            let h = 0.01;
            let b = SpatialBushing {
                rest_translation_a: if torsion { Z } else { Vec3::new(1., 0., 0.) },
                translation_stiffness: if torsion { Z } else { Vec3::new(k, 0., 0.) },
                translation_damping: if torsion {
                    Z
                } else {
                    Vec3::new(damping, 0., 0.)
                },
                rotation_stiffness: if torsion { Vec3::new(0., 0., k) } else { Z },
                rotation_damping: if torsion {
                    Vec3::new(0., 0., damping)
                } else {
                    Z
                },
                ..element(frame(0, Z, Z), frame(1, Z, Z))
            };
            let mut state = FreeFloatingMotionState {
                base: FreeFloatingBaseState::stationary(Se3::identity()),
                joint_position: vec![0.1],
                joint_velocity: vec![-0.4],
            };
            let mut q = 0.1;
            let mut v = -0.4;
            let mut maximum_evaluations = 0;
            for _ in 0..40 {
                let mid = (v - 0.5 * h * k * q / mu)
                    / (1. + 0.5 * h * damping / mu + 0.25 * h * h * k / mu);
                q += h * mid;
                v = 2. * mid - v;
                let result = step_free_floating_bushings_implicit(
                    &model,
                    &state,
                    held(&[0.], &[ZERO; 2]),
                    &[b],
                    h,
                    ImplicitSpatialConfig {
                        residual_tolerance: 1e-8,
                        ..Default::default()
                    },
                    || false,
                )
                .unwrap();
                maximum_evaluations = maximum_evaluations.max(result.residual_evaluations);
                let step = result.spatial;
                close(step.motion.state_after.joint_position[0], q, 2e-10);
                close(step.motion.state_after.joint_velocity[0], v, 3e-8);
                close(step.dissipated_energy, h * damping * mid * mid, 3e-7);
                close(step.stored_energy_after, 0.5 * k * q * q, 3e-5);
                assert!(
                    step.mechanical_work_residual.abs() < 3e-7,
                    "{}",
                    step.mechanical_work_residual
                );
                balance(&step.motion, 1e-10);
                state = step.motion.state_after;
            }
            eprintln!(
                "stiff bushing torsion={torsion} damping={damping}: max evaluations {maximum_evaluations}"
            );
        }
    }
}

fn state_error(a: &FreeFloatingMotionState, b: &FreeFloatingMotionState) -> f64 {
    let mut sum = (a.base.world_from_base.translation() - b.base.world_from_base.translation())
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
        sum += (x - y).powi(2);
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
        sum += (x - y).powi(2);
    }
    sum.sqrt()
}

#[test]
fn mixed_compliant_tree_converges_with_original_external_impulses_and_total_storage() {
    let (model, start, bs) = mixed();
    let springs = [AxialConnection {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::new(0.1, 0.2, 0.),
        },
        b: MaterialPoint {
            link: 3,
            point_body: Vec3::new(-0.1, 0.1, 0.2),
        },
        rest_length_m: 0.8,
        stiffness_n_per_m: 8.,
        damping_ns_per_m: 0.4,
    }];
    let points = [WorldPointLoad {
        point: MaterialPoint {
            link: 2,
            point_body: Vec3::new(0.3, -0.2, 0.1),
        },
        force_world: Vec3::new(1., -0.5, 0.3),
        couple_world: Vec3::new(0.1, 0.2, -0.1),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: Vec3::new(0., -0.2, 0.1),
            ..held(&[0.1, -0.2, 0.05], &[ZERO; 4]).held
        },
        world_points: &points,
        connections: &springs,
    };
    let duration = 0.12;
    let mut reference = start.clone();
    for _ in 0..2048 {
        reference =
            step_free_floating_bushings(&model, &reference, loads, &bs, duration / 2048., || false)
                .unwrap()
                .motion
                .state_after;
    }
    let mut errors = Vec::new();
    let mut defects = Vec::new();
    for count in [8, 16, 32] {
        let mut state = start.clone();
        let mut defect = 0.;
        for _ in 0..count {
            let step = step_free_floating_bushings_implicit(
                &model,
                &state,
                loads,
                &bs,
                duration / f64::from(count),
                ImplicitSpatialConfig::default(),
                || false,
            )
            .unwrap()
            .spatial;
            balance(&step.motion, 8e-12);
            assert!(step.dissipated_energy >= 0.);
            let axial = |s: &FreeFloatingMotionState| {
                let fk = forward_kinematics(
                    &model,
                    BaseState::stationary(Se3::identity()),
                    &s.joint_position,
                    &s.joint_velocity,
                )
                .unwrap();
                let a = fk.world_from_link[1]
                    .transform_point(springs[0].a.point_body)
                    .unwrap();
                let b = fk.world_from_link[3]
                    .transform_point(springs[0].b.point_body)
                    .unwrap();
                0.5 * springs[0].stiffness_n_per_m
                    * ((b - a).norm() - springs[0].rest_length_m).powi(2)
            };
            close(
                step.stored_energy_before,
                reference_energy(&model, &state, &bs) + axial(&state),
                1e-12,
            );
            close(
                step.stored_energy_after,
                reference_energy(&model, &step.motion.state_after, &bs)
                    + axial(&step.motion.state_after),
                1e-12,
            );
            defect += step.mechanical_work_residual;
            state = step.motion.state_after;
        }
        errors.push(state_error(&state, &reference));
        defects.push(defect.abs());
    }
    eprintln!("mixed bushing tree errors {errors:?}; work defects {defects:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8);
    }
    for pair in defects.windows(2) {
        assert!(pair[0] / pair[1] > 3.5);
    }
}

#[test]
fn empty_bushings_reuse_both_original_integrators_and_zero_axes_allow_unrestricted_rotation() {
    let (model, state, bs) = mixed();
    let loads = held(&[0.; 3], &[ZERO; 4]);
    assert_eq!(
        step_free_floating_bushings(&model, &state, loads, &[], 0.002, || false).unwrap(),
        step_free_floating_spatial(&model, &state, loads, 0.002, || false).unwrap()
    );
    assert_eq!(
        step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &[],
            0.002,
            Default::default(),
            || false
        )
        .unwrap(),
        step_free_floating_spatial_implicit(
            &model,
            &state,
            loads,
            0.002,
            Default::default(),
            || false
        )
        .unwrap()
    );
    let b = SpatialBushing {
        b: frame(2, Vec3::new(0., 0., core::f64::consts::PI), Z),
        ..element(bs[0].a, bs[0].b)
    };
    let response = evaluate_bushing_forces(&model, &state, &[b], || false).unwrap();
    assert_eq!(response.body_wrench, vec![ZERO; 4]);
    assert_eq!(response.stored_energy, 0.);
    assert_eq!(response.internal_power, 0.);
}

#[test]
fn invalid_material_chart_budget_and_late_cancellation_refuse_without_mutation() {
    let (model, state, bs) = mixed();
    let before = state.clone();
    let loads = held(&[0.; 3], &[ZERO; 4]);
    let invalid = SpatialBushing {
        rotation_damping: Vec3::new(0., -1., 0.),
        ..bs[0]
    };
    assert!(evaluate_bushing_forces(&model, &state, &[invalid], || false).is_err());
    let bad_link = SpatialBushing {
        a: MaterialFrame {
            link: 99,
            ..bs[0].a
        },
        ..bs[0]
    };
    assert!(evaluate_bushing_forces(&model, &state, &[bad_link], || false).is_err());
    let fk = forward_kinematics(
        &model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let frame_a = fk.world_from_link[bs[0].a.link]
        .rotation()
        .compose(bs[0].a.link_from_attachment.rotation())
        .unwrap();
    let chart = SpatialBushing {
        b: MaterialFrame {
            link: bs[0].b.link,
            link_from_attachment: Se3::from_parts(
                fk.world_from_link[bs[0].b.link]
                    .rotation()
                    .inverse()
                    .compose(frame_a)
                    .unwrap()
                    .compose(rot(Vec3::new(
                        0.,
                        0.,
                        MAX_BUSHING_ROTATION_RADIANS + 0.5e-6,
                    )))
                    .unwrap(),
                Z,
            )
            .unwrap(),
        },
        ..bs[0]
    };
    assert!(matches!(
        evaluate_bushing_forces(&model, &state, &[chart], || false),
        Err(FreeFloatingMotionError::InvalidLoad(_))
    ));
    assert!(matches!(
        step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &bs,
            0.02,
            ImplicitSpatialConfig {
                maximum_evaluations: 2,
                ..Default::default()
            },
            || false
        ),
        Err(ImplicitSpatialError::EvaluationBudget { .. })
    ));
    let mut polls = 0;
    let accepted = step_free_floating_bushings_implicit(
        &model,
        &state,
        loads,
        &bs,
        0.002,
        Default::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    let cancelled = step_free_floating_bushings_implicit(
        &model,
        &state,
        loads,
        &bs,
        0.002,
        Default::default(),
        || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(
        cancelled,
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Cancelled
        ))
    ));
    assert_eq!(state, before);
    assert_eq!(
        accepted,
        step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &bs,
            0.002,
            Default::default(),
            || false
        )
        .unwrap()
    );
}
