//! Physical spherical closure, mixed loops, native recoil and reaction work.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, BilateralLoopStep, DistanceLoop,
    PointCoincidenceLoop, step_free_floating_bilateral_loops,
    step_free_floating_bilateral_loops_with_bushings,
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
    apply_free_floating_impulses, forward_kinematics, free_floating_momentum,
};

const Z: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};
const ZERO: Wrench = Wrench::new(Z, Z);
const AXES: [Vec3; 3] = [
    Vec3 {
        x: 1.,
        y: 0.,
        z: 0.,
    },
    Vec3 {
        x: 0.,
        y: 1.,
        z: 0.,
    },
    Vec3 {
        x: 0.,
        y: 0.,
        z: 1.,
    },
];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; error {:.4e}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!(
        (a - b).norm() <= tolerance,
        "{a:?} != {b:?}; error {}",
        (a - b).norm()
    );
}
fn rot(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn pose(v: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), v).unwrap()
}
fn point(link: usize, p: Vec3) -> MaterialPoint {
    MaterialPoint {
        link,
        point_body: p,
    }
}
fn ball(a: MaterialPoint, b: MaterialPoint) -> BilateralLoop {
    BilateralLoop::PointCoincidence(PointCoincidenceLoop { a, b })
}
fn inertia(mass: f64, center: Vec3) -> SpatialInertia {
    let mut tensor = [0.3, 0.025, -0.012, 0.025, 0.4, 0.018, -0.012, 0.018, 0.5];
    for value in &mut tensor {
        *value *= mass;
    }
    SpatialInertia::new(mass, center, Mat3 { m: tensor }).unwrap()
}
fn state(q: Vec<f64>, v: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: q,
        joint_velocity: v,
    }
}
fn loads<'a>(effort: &'a [f64], wrench: &'a [Wrench]) -> SpatialMotionLoads<'a> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: effort,
            external_body_wrench: wrench,
            gravity_world: Z,
        },
        connections: &[],
        world_points: &[],
    }
}
fn account(s: &BilateralLoopStep, tolerance: f64) {
    assert!(s.position_residual_after_m.iter().all(|x| x.abs() <= 1e-9));
    assert!(
        s.velocity_residual_after_mps
            .iter()
            .all(|x| x.abs() <= 1e-9)
    );
    for p in &s.point_residual_after {
        assert!(p.separation_world_m.norm() <= 1e-9);
        assert!(p.relative_velocity_world_mps.norm() <= 1e-9);
    }
    vector(s.linear_momentum_residual_world, Z, tolerance);
    vector(s.angular_momentum_residual_world, Z, tolerance);
    vector(s.reaction_linear_impulse_world, Z, tolerance);
    vector(s.reaction_angular_impulse_world, Z, tolerance);
    vector(s.linear_momentum_balance_residual_world, Z, tolerance);
    vector(s.angular_momentum_balance_residual_world, Z, tolerance);
    close(
        s.mechanical_work_residual,
        s.constraint_impulse_work
            + s.midpoint.mechanical_work_residual
            + s.initial_reaction.kinetic_work_residual
            + s.final_reaction.kinetic_work_residual,
        tolerance,
    );
}

// Independent parallel-axis sum, with no native articulated/composite-inertia map.
fn composite(parts: &[(SpatialInertia, Vec3)]) -> SpatialInertia {
    let mass: f64 = parts.iter().map(|(i, _)| i.mass()).sum();
    let mut center = Z;
    for (i, x) in parts {
        center = center + (*x + i.center_of_mass()).scale(i.mass() / mass);
    }
    let mut tensor = [0.; 9];
    for (i, x) in parts {
        let d = *x + i.center_of_mass() - center;
        let a = [d.x, d.y, d.z];
        for row in 0..3 {
            for col in 0..3 {
                tensor[3 * row + col] += i.inertia_com().m[3 * row + col]
                    + i.mass() * (if row == col { d.dot(d) } else { 0. } - a[row] * a[col]);
            }
        }
    }
    SpatialInertia::new(mass, center, Mat3 { m: tensor }).unwrap()
}

fn ball_reduction() -> (
    ArticulatedModel,
    ArticulatedModel,
    FreeFloatingMotionState,
    BilateralLoop,
) {
    let pivot = Vec3::new(0.35, -0.25, 0.2);
    let root_inertia = inertia(3., Vec3::new(0.1, -0.15, 0.05));
    let slides = [
        inertia(0.4, Vec3::new(0.04, -0.03, 0.02)),
        inertia(0.5, Vec3::new(-0.02, 0.03, 0.01)),
        inertia(0.6, Vec3::new(0.01, 0.02, -0.04)),
    ];
    let rotating = [
        inertia(0.7, Vec3::new(0.12, 0.04, -0.03)),
        inertia(0.8, Vec3::new(-0.03, 0.16, 0.06)),
        inertia(1.2, Vec3::new(0.3, -0.2, 0.1)),
    ];
    let mut links = vec![Link::new(
        "base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        root_inertia,
    )];
    for i in 0..3 {
        links.push(Link::new(
            format!("slide_{i}"),
            Some(i),
            if i == 0 { pose(pivot) } else { Se3::identity() },
            JointModel::prismatic(AXES[i], None).unwrap(),
            slides[i],
        ));
    }
    for i in 0..3 {
        links.push(Link::new(
            format!("rotation_{i}"),
            Some(i + 3),
            Se3::identity(),
            JointModel::revolute(AXES[i], None).unwrap(),
            rotating[i],
        ));
    }
    let combined = composite(&[
        (root_inertia, Z),
        (slides[0], pivot),
        (slides[1], pivot),
        (slides[2], pivot),
    ]);
    let mut reduced = vec![Link::new(
        "combined_base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        combined,
    )];
    for i in 0..3 {
        reduced.push(Link::new(
            format!("rotation_{i}"),
            Some(i),
            if i == 0 { pose(pivot) } else { Se3::identity() },
            JointModel::revolute(AXES[i], None).unwrap(),
            rotating[i],
        ));
    }
    let mut s = state(
        vec![0., 0., 0., 0.3, -0.4, 0.2],
        vec![0., 0., 0., 0.5, -0.3, 0.6],
    );
    s.base = FreeFloatingBaseState::new(
        Se3::from_parts(rot(Vec3::new(0.2, -0.1, 0.3)), Vec3::new(0.1, -0.2, 0.4)).unwrap(),
        Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.5, -0.1, 0.2)),
    );
    (
        ArticulatedModel::new(links).unwrap(),
        ArticulatedModel::new(reduced).unwrap(),
        s,
        ball(point(0, pivot), point(6, Z)),
    )
}

fn state_error(a: &FreeFloatingMotionState, b: &FreeFloatingMotionState, skip: usize) -> f64 {
    let mut e = (a.base.world_from_base.translation() - b.base.world_from_base.translation())
        .norm()
        .powi(2)
        + (a.base.twist_body.linear - b.base.twist_body.linear)
            .norm()
            .powi(2)
        + (a.base.twist_body.angular - b.base.twist_body.angular)
            .norm()
            .powi(2);
    for (x, y) in a
        .base
        .world_from_base
        .rotation()
        .matrix()
        .m
        .iter()
        .zip(b.base.world_from_base.rotation().matrix().m)
    {
        e += (x - y).powi(2);
    }
    for (x, y) in a.joint_position[skip..].iter().zip(&b.joint_position) {
        e += (x - y).powi(2);
    }
    for (x, y) in a.joint_velocity[skip..].iter().zip(&b.joint_velocity) {
        e += (x - y).powi(2);
    }
    e.sqrt()
}

#[test]
fn spherical_closure_converges_to_an_independently_reduced_three_rotation_mechanism() {
    let (model, reduced, initial, constraint) = ball_reduction();
    let mut reference = FreeFloatingMotionState {
        base: initial.base,
        joint_position: initial.joint_position[3..].to_vec(),
        joint_velocity: initial.joint_velocity[3..].to_vec(),
    };
    let original = free_floating_momentum(
        &model,
        initial.base,
        &initial.joint_position,
        &initial.joint_velocity,
    )
    .unwrap();
    let reduced_initial = free_floating_momentum(
        &reduced,
        reference.base,
        &reference.joint_position,
        &reference.joint_velocity,
    )
    .unwrap();
    vector(
        original.linear_momentum_world,
        reduced_initial.linear_momentum_world,
        2e-13,
    );
    vector(
        original.angular_momentum_world,
        reduced_initial.angular_momentum_world,
        2e-13,
    );
    close(
        original.kinetic_energy,
        reduced_initial.kinetic_energy,
        2e-13,
    );
    let force = WorldPointLoad {
        point: point(6, Vec3::new(0.12, 0.07, -0.04)),
        force_world: Vec3::new(0.5, -0.3, 0.2),
        couple_world: Vec3::new(0.2, 0.1, -0.25),
    };
    let reference_force = WorldPointLoad {
        point: point(3, force.point.point_body),
        ..force
    };
    let efforts = [0., 0., 0., 0.3, -0.2, 0.4];
    let reference_efforts = [0.3, -0.2, 0.4];
    let body = [ZERO; 7];
    let reference_body = [ZERO; 4];
    let forces = [force];
    let reference_forces = [reference_force];
    let mut load = loads(&efforts, &body);
    load.held.gravity_world = Vec3::new(-0.1, -0.4, -0.2);
    load.world_points = &forces;
    let mut reference_load = loads(&reference_efforts, &reference_body);
    reference_load.held.gravity_world = load.held.gravity_world;
    reference_load.world_points = &reference_forces;
    let duration = 0.24;
    for _ in 0..512 {
        reference = step_free_floating_spatial_implicit(
            &reduced,
            &reference,
            reference_load,
            duration / 512.,
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
    for count in [6, 12, 24] {
        let mut s = initial.clone();
        let mut defect = 0.;
        for _ in 0..count {
            let step = step_free_floating_bilateral_loops(
                &model,
                &s,
                load,
                &[constraint],
                duration / f64::from(count),
                Default::default(),
                || false,
            )
            .unwrap();
            account(&step, 4e-11);
            for &q in &step.state_after.joint_position[..3] {
                close(q, 0., 1e-10);
            }
            for &v in &step.state_after.joint_velocity[..3] {
                close(v, 0., 1e-10);
            }
            defect += step.mechanical_work_residual;
            s = step.state_after;
        }
        errors.push(state_error(&s, &reference, 3));
        defects.push(defect.abs());
        assert!((s.joint_position[5] - initial.joint_position[5]).abs() > 0.1);
    }
    eprintln!("ball-joint independent reduced errors {errors:?}; mechanical defects {defects:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.4, "{errors:?}");
    }
}

fn crank() -> (ArticulatedModel, FreeFloatingMotionState, BilateralLoop) {
    let mut links = vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3., Z),
        ),
        Link::new(
            "crank",
            Some(0),
            Se3::identity(),
            JointModel::revolute(AXES[2], None).unwrap(),
            inertia(1., Vec3::new(0.5, 0., 0.)),
        ),
    ];
    for i in 0..3 {
        links.push(Link::new(
            format!("slide_{i}"),
            Some(if i == 0 { 0 } else { i + 1 }),
            Se3::identity(),
            JointModel::prismatic(AXES[i], None).unwrap(),
            inertia(0.6, Z),
        ));
    }
    let q: f64 = 0.4;
    let v = 0.7;
    let mut s = state(
        vec![q, q.cos(), q.sin(), 0.],
        vec![v, -q.sin() * v, q.cos() * v, 0.],
    );
    s.base.twist_body = Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.4, -0.2, 0.1));
    (
        ArticulatedModel::new(links).unwrap(),
        s,
        ball(point(1, AXES[0]), point(4, Z)),
    )
}

#[test]
fn rotating_material_point_closure_and_finite_bushing_share_the_nonlinear_stage() {
    let (model, initial, constraint) = crank();
    let element = SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 1,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Z,
        translation_stiffness: Z,
        translation_damping: Z,
        rotation_stiffness: Vec3::new(1., 2., 3.),
        rotation_damping: Vec3::new(0.2, 0.3, 0.5),
    };
    let mut s = initial.clone();
    let effort = [0.3, 0., 0., 0.];
    let body = [ZERO; 5];
    let force = [WorldPointLoad {
        point: point(4, Z),
        force_world: Vec3::new(0.3, -0.2, 0.1),
        couple_world: Z,
    }];
    let load = SpatialMotionLoads {
        world_points: &force,
        ..loads(&effort, &body)
    };
    let plain = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &[constraint],
        0.01,
        Default::default(),
        || false,
    )
    .unwrap();
    let empty = step_free_floating_bilateral_loops_with_bushings(
        &model,
        &s,
        load,
        &[constraint],
        &[],
        0.01,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(plain, empty);
    let mut loss = 0.;
    let mut work = 0.;
    let mut maximum_gap = 0.0_f64;
    for _ in 0..30 {
        let result = step_free_floating_bilateral_loops_with_bushings(
            &model,
            &s,
            load,
            &[constraint],
            &[element],
            0.01,
            Default::default(),
            || false,
        )
        .unwrap();
        account(&result, 2e-11);
        close(
            result.midpoint.stored_energy_before,
            1.5 * s.joint_position[0].powi(2),
            2e-13,
        );
        let q = result.state_after.joint_position[0];
        let v = result.state_after.joint_velocity[0];
        close(result.state_after.joint_position[1], q.cos(), 1e-10);
        close(result.state_after.joint_position[2], q.sin(), 1e-10);
        close(result.state_after.joint_velocity[1], -q.sin() * v, 1e-10);
        close(result.state_after.joint_velocity[2], q.cos() * v, 1e-10);
        let actual = forward_kinematics(
            &model,
            BaseState::prescribed(
                result.state_after.base.world_from_base,
                result.state_after.base.twist_body,
                Twist::zero(),
            ),
            &result.state_after.joint_position,
            &result.state_after.joint_velocity,
        )
        .unwrap();
        let gap = actual.world_from_link[4].translation()
            - actual.world_from_link[1].transform_point(AXES[0]).unwrap();
        vector(
            gap,
            result.point_residual_after[0].separation_world_m,
            4e-16,
        );
        maximum_gap = maximum_gap.max(gap.norm());
        loss += result.midpoint.dissipated_energy;
        work += result.mechanical_work_residual;
        s = result.state_after;
    }
    assert!(loss > 0.01);
    assert!((s.joint_position[0] - initial.joint_position[0]).abs() > 0.1);
    eprintln!(
        "nonlinear material-frame ball loop: max actual gap={maximum_gap:e}, bushing loss={loss:e}, mechanical defect={work:e}"
    );
}

#[test]
fn two_ball_loops_and_a_distance_share_one_full_mobility_and_force_balance() {
    let mut links = vec![Link::new(
        "base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(3., Z),
    )];
    for branch in 0..2 {
        let first = links.len();
        for (i, axis) in AXES.into_iter().enumerate() {
            links.push(Link::new(
                format!("branch_{branch}_{i}"),
                Some(if i == 0 { 0 } else { first + i - 1 }),
                Se3::identity(),
                JointModel::prismatic(axis, None).unwrap(),
                inertia(1., Z),
            ));
        }
    }
    links.push(Link::new(
        "distance_slider",
        Some(0),
        pose(AXES[2]),
        JointModel::prismatic(AXES[2], None).unwrap(),
        inertia(2., Z),
    ));
    let model = ArticulatedModel::new(links).unwrap();
    let s = state(vec![0.; 7], vec![0.; 7]);
    let constraints = [
        ball(point(0, Z), point(3, Z)),
        BilateralLoop::Distance(DistanceLoop {
            a: point(0, Z),
            b: point(7, Z),
            length_m: 1.,
        }),
        ball(point(0, Z), point(6, Z)),
    ];
    let effort = [0.; 7];
    let body = [ZERO; 8];
    let force = 2.7;
    let h = 0.03;
    let applied = [WorldPointLoad {
        point: point(6, Z),
        force_world: AXES[2].scale(force),
        couple_world: Z,
    }];
    let load = SpatialMotionLoads {
        world_points: &applied,
        ..loads(&effort, &body)
    };
    let result = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &constraints,
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&result, 2e-12);
    assert_eq!(result.rows.len(), 7);
    assert_eq!(result.rows[3].component, None);
    assert_eq!(result.rows[4].constraint, 2);
    vector(
        result.state_after.base.twist_body.linear,
        AXES[2].scale(force * h / 11.),
        2e-13,
    );
    vector(
        result.state_after.base.world_from_base.translation(),
        AXES[2].scale(0.5 * force * h * h / 11.),
        2e-13,
    );
    vector(result.state_after.base.twist_body.angular, Z, 2e-13);
    for &v in &result.state_after.joint_velocity {
        close(v, 0., 2e-12);
    }
    let reordered = [constraints[2], constraints[0], constraints[1]];
    let permuted = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &reordered,
        h,
        Default::default(),
        || false,
    )
    .unwrap();
    close(
        state_error(&result.state_after, &permuted.state_after, 0),
        0.,
        2e-12,
    );
    close(result.constraint_impulse_work, 0., 2e-13);
}

#[test]
fn admitted_finite_gap_reaction_includes_the_true_geometric_couple_and_its_work() {
    let (model, _, mut s, constraint) = ball_reduction();
    s.joint_position[..3].copy_from_slice(&[2e-6, -1e-6, 0.5e-6]);
    let effort = [0., 0., 0., 0.3, -0.2, 0.4];
    let body = [ZERO; 7];
    let applied = [WorldPointLoad {
        point: point(6, Vec3::new(0.2, -0.1, 0.15)),
        force_world: Vec3::new(1.3, -0.7, 0.9),
        couple_world: Vec3::new(0.2, 0.1, -0.3),
    }];
    let load = SpatialMotionLoads {
        world_points: &applied,
        ..loads(&effort, &body)
    };
    let config = BilateralLoopConfig {
        position_tolerance_m: 1e-5,
        velocity_tolerance_mps: 1e-5,
        ..Default::default()
    };
    let result =
        step_free_floating_bilateral_loops(&model, &s, load, &[constraint], 0.02, config, || false)
            .unwrap();
    account(&result, 5e-11);
    let BilateralLoop::PointCoincidence(c) = constraint else {
        unreachable!()
    };
    let before = forward_kinematics(
        &model,
        BaseState::prescribed(s.base.world_from_base, s.base.twist_body, Twist::zero()),
        &s.joint_position,
        &s.joint_velocity,
    )
    .unwrap();
    let a = before.world_from_link[c.a.link];
    let b = before.world_from_link[c.b.link];
    let lambda = Vec3::new(
        result.initial_impulse_ns[0],
        result.initial_impulse_ns[1],
        result.initial_impulse_ns[2],
    );
    let force_b_world = a.rotation().rotate(lambda).unwrap();
    let couple_world = result.point_residual_before[0]
        .separation_world_m
        .cross(force_b_world)
        .scale(-1.);
    assert!(couple_world.norm() > 1e-9, "{couple_world:?}");
    vector(
        result.initial_geometric_couple_impulse_world[0],
        couple_world,
        1e-20,
    );
    let fa = a
        .rotation()
        .inverse()
        .rotate(force_b_world.scale(-1.))
        .unwrap();
    let fb = b.rotation().inverse().rotate(force_b_world).unwrap();
    let mut expected = vec![ZERO; model.link_count()];
    expected[c.a.link] = Wrench::new(
        c.a.point_body.cross(fa) + a.rotation().inverse().rotate(couple_world).unwrap(),
        fa,
    );
    expected[c.b.link] = Wrench::new(c.b.point_body.cross(fb), fb);
    let independent = apply_free_floating_impulses(
        &model,
        s.base,
        &s.joint_position,
        &s.joint_velocity,
        &[0.; 6],
        &expected,
    )
    .unwrap();
    vector(
        independent.base_after.twist_body.angular,
        result.initial_reaction.base_after.twist_body.angular,
        3e-14,
    );
    vector(
        independent.base_after.twist_body.linear,
        result.initial_reaction.base_after.twist_body.linear,
        3e-14,
    );
    for (x, y) in independent
        .joint_velocity_after
        .iter()
        .zip(&result.initial_reaction.joint_velocity_after)
    {
        close(*x, *y, 3e-14);
    }
    close(
        independent.impulse_work,
        result.initial_reaction.impulse_work,
        2e-14,
    );
    let midpoint_omega = a
        .rotation()
        .rotate(
            (s.base.twist_body.angular + result.initial_reaction.base_after.twist_body.angular)
                .scale(0.5),
        )
        .unwrap();
    close(
        result.initial_geometric_couple_work,
        couple_world.dot(midpoint_omega),
        1e-20,
    );
    expected[c.a.link].torque = c.a.point_body.cross(fa);
    let missing_couple = apply_free_floating_impulses(
        &model,
        s.base,
        &s.joint_position,
        &s.joint_velocity,
        &[0.; 6],
        &expected,
    )
    .unwrap();
    vector(
        missing_couple.external_angular_impulse_world,
        couple_world.scale(-1.),
        2e-16,
    );
    assert!(missing_couple.external_angular_impulse_world.norm() > 1e-9);
    vector(
        result.initial_reaction.external_angular_impulse_world,
        Z,
        2e-16,
    );
    eprintln!(
        "admitted-gap geometric couple={:e} N m s; its retained impulse work={:e} J",
        couple_world.norm(),
        result.initial_geometric_couple_work
    );
}

#[test]
fn point_loop_motion_and_material_reactions_are_scene_covariant() {
    let (model, s, constraint) = crank();
    let effort = [0.3, 0., 0., 0.];
    let body = [ZERO; 5];
    let applied = [WorldPointLoad {
        point: point(4, Z),
        force_world: Vec3::new(0.7, -0.4, 0.2),
        couple_world: Vec3::new(0.1, 0.2, -0.3),
    }];
    let mut load = SpatialMotionLoads {
        world_points: &applied,
        ..loads(&effort, &body)
    };
    load.held.gravity_world = Vec3::new(0.1, -0.3, 0.2);
    let scene = Se3::from_parts(rot(Vec3::new(0.5, -0.3, 0.7)), Vec3::new(2., -3., 4.)).unwrap();
    let mut transformed = s.clone();
    transformed.base.world_from_base = scene.compose(s.base.world_from_base).unwrap();
    let rotated_loads = [WorldPointLoad {
        force_world: scene.rotation().rotate(applied[0].force_world).unwrap(),
        couple_world: scene.rotation().rotate(applied[0].couple_world).unwrap(),
        ..applied[0]
    }];
    let mut transformed_load = SpatialMotionLoads {
        world_points: &rotated_loads,
        ..load
    };
    transformed_load.held.gravity_world = scene.rotation().rotate(load.held.gravity_world).unwrap();
    let a = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &[constraint],
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    let b = step_free_floating_bilateral_loops(
        &model,
        &transformed,
        transformed_load,
        &[constraint],
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    account(&a, 2e-11);
    account(&b, 2e-11);
    let mut expected = a.state_after.clone();
    expected.base.world_from_base = scene.compose(expected.base.world_from_base).unwrap();
    close(state_error(&expected, &b.state_after, 0), 0., 3e-11);
    for (x, y) in a.initial_impulse_ns.iter().zip(&b.initial_impulse_ns) {
        close(*x, *y, 3e-11);
    }
    for (x, y) in a.final_impulse_ns.iter().zip(&b.final_impulse_ns) {
        close(*x, *y, 3e-11);
    }
    close(a.constraint_impulse_work, b.constraint_impulse_work, 3e-12);
    close(
        a.mechanical_work_residual,
        b.mechanical_work_residual,
        3e-12,
    );
}

#[test]
fn point_norm_physical_tangency_and_mixed_rank_are_admitted_honestly() {
    let (model, _, initial, constraint) = ball_reduction();
    let effort = [0.; 6];
    let body = [ZERO; 7];
    let load = loads(&effort, &body);
    let run = |s: &FreeFloatingMotionState, constraints: &[BilateralLoop]| {
        step_free_floating_bilateral_loops(
            &model,
            s,
            load,
            constraints,
            0.01,
            Default::default(),
            || false,
        )
    };
    let mut separated = initial.clone();
    separated.joint_position[..3].copy_from_slice(&[0.8e-9, 0.8e-9, 0.]);
    assert!(matches!(
        run(&separated, &[constraint]),
        Err(BilateralLoopError::Position {
            stage: "initial",
            ..
        })
    ));
    let mut spinning_gap = initial.clone();
    spinning_gap.joint_position[0] = 0.9e-9;
    spinning_gap.base.twist_body.angular = Vec3::new(0., 0., 10.);
    assert!(matches!(
        run(&spinning_gap, &[constraint]),
        Err(BilateralLoopError::PointVelocity {
            stage: "initial",
            ..
        })
    ));
    let mut nontangent = initial.clone();
    nontangent.joint_velocity[1] = 1e-3;
    assert!(matches!(
        run(&nontangent, &[constraint]),
        Err(BilateralLoopError::Velocity {
            stage: "initial",
            ..
        })
    ));
    assert!(matches!(
        run(&initial, &[constraint, constraint]),
        Err(BilateralLoopError::Mobility {
            stage: "initial",
            ..
        })
    ));
    let pivot = Vec3::new(0.35, -0.25, 0.2);
    let redundant_distance = BilateralLoop::Distance(DistanceLoop {
        a: point(0, pivot),
        b: point(6, AXES[0]),
        length_m: 1.,
    });
    assert!(matches!(
        run(&initial, &[constraint, redundant_distance]),
        Err(BilateralLoopError::Mobility {
            stage: "initial",
            constraint: 3,
            ..
        })
    ));
    let imposed = ball(point(3, Z), point(6, Z));
    assert!(matches!(
        run(&initial, &[imposed]),
        Err(BilateralLoopError::Mobility {
            stage: "initial",
            ..
        })
    ));
    assert!(matches!(
        run(&initial, &[constraint; 11]),
        Err(BilateralLoopError::Invalid(_))
    ));
    assert!(matches!(
        run(&initial, &[]),
        Err(BilateralLoopError::Invalid(_))
    ));
    assert!(matches!(
        run(&initial, &[ball(point(0, Z), point(0, AXES[0]))]),
        Err(BilateralLoopError::Invalid(_))
    ));
    assert!(matches!(
        run(
            &initial,
            &[ball(point(0, Z), point(6, Vec3::new(f64::NAN, 0., 0.)))]
        ),
        Err(_)
    ));
}

#[test]
fn point_loop_cancellation_final_budget_and_deterministic_retry_are_atomic() {
    let (model, s, constraint) = crank();
    let effort = [0.3, 0., 0., 0.];
    let body = [ZERO; 5];
    let load = loads(&effort, &body);
    let untouched = s.clone();
    let mut polls = 0usize;
    let baseline = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &[constraint],
        0.02,
        Default::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        let error = step_free_floating_bilateral_loops(
            &model,
            &s,
            load,
            &[constraint],
            0.02,
            Default::default(),
            || {
                calls += 1;
                calls >= stop
            },
        )
        .unwrap_err();
        assert!(
            matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Cancelled
                ))
            ),
            "{error:?}"
        );
        assert_eq!(s, untouched);
    }
    let mut config = BilateralLoopConfig::default();
    config.nonlinear.maximum_evaluations = baseline.residual_evaluations - 1;
    assert!(matches!(
        step_free_floating_bilateral_loops(&model, &s, load, &[constraint], 0.02, config, || false),
        Err(BilateralLoopError::Motion(
            ImplicitSpatialError::EvaluationBudget { .. }
        ))
    ));
    let retried = step_free_floating_bilateral_loops(
        &model,
        &s,
        load,
        &[constraint],
        0.02,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(baseline, retried);
    assert_eq!(s, untouched);
    eprintln!(
        "point-loop bounded solve evaluations={}, cancellation polls={polls}",
        baseline.residual_evaluations
    );
}
