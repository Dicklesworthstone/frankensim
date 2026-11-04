//! G0/G1/G3/G4: physical closed mechanisms and independently driven supports.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::prescribed::motion::{
    PrescribedBilateralLoopStep, step_with_bushings, step_with_bushings_and_forcing,
    step_with_coupled_bushings_and_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, DistanceLoop, PointCoincidenceLoop,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitPrescribedSpatialError, step_prescribed_base_bushings_implicit_with_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_kinematics,
};

const T0: f64 = 4.0;
const A: f64 = 0.23;
const ALPHA: f64 = -0.17;
const Q0: [f64; 2] = [0.3, -0.2];
const V0: [f64; 2] = [-0.1, 0.4];
const HELD: [f64; 5] = [0.3, -0.2, 0.0, 0.0, 0.0];
const WRENCHES: [Wrench; 6] = [Wrench::new(
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
); 6];
const K: [[f64; 2]; 2] = [[12.0, 4.0], [4.0, 8.0]];
const C: [[f64; 2]; 2] = [[0.6, 0.1], [0.1, 0.4]];
const DRAG: [f64; 2] = [70.0, 45.0];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:.16e} != {b:.16e}; tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() < tolerance, "{a:?} != {b:?}");
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn pose(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0., 0.)).unwrap()
}
fn inertia(mass: f64, spin: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [spin, 0., 0., 0., spin, 0., 0., 0., spin],
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
fn config() -> BilateralLoopConfig {
    BilateralLoopConfig {
        nonlinear: ImplicitSpatialConfig {
            residual_tolerance: 1e-12,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn support(time: f64) -> BaseState {
    let t = time - T0;
    BaseState::prescribed(
        Se3::from_parts(
            rotation(Vec3::new(0.2 + 0.3 * t + 0.5 * ALPHA * t * t, 0., 0.)),
            Vec3::new(2.0 + 0.4 * t + 0.5 * A * t * t, -0.7, 0.4),
        )
        .unwrap(),
        Twist::new(
            Vec3::new(0.3 + ALPHA * t, 0., 0.),
            Vec3::new(0.4 + A * t, 0., 0.),
        ),
        Twist::new(Vec3::new(ALPHA, 0., 0.), Vec3::new(A, 0., 0.)),
    )
}
fn path(time: f64, h: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: support(time + 0.5 * h),
        endpoint: support(time + h),
    }
}
fn model(ball: bool) -> ArticulatedModel {
    let slide = |axis| JointModel::prismatic(axis, None).unwrap();
    let mut links = vec![
        Link::new(
            "carrier",
            None,
            pose(0.),
            JointModel::FIXED,
            inertia(3., 3.),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(1.),
            slide(Vec3::new(1., 0., 0.)),
            inertia(1., 2.),
        ),
        Link::new(
            "rotor",
            Some(1),
            pose(0.),
            JointModel::revolute(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(1., 2.),
        ),
        Link::new(
            "follower x",
            Some(0),
            pose(2.),
            slide(Vec3::new(1., 0., 0.)),
            inertia(1., 3.),
        ),
    ];
    if ball {
        links.push(Link::new(
            "follower y",
            Some(3),
            pose(0.),
            slide(Vec3::new(0., 1., 0.)),
            inertia(1., 3.),
        ));
        links.push(Link::new(
            "follower z",
            Some(4),
            pose(0.),
            slide(Vec3::new(0., 0., 1.)),
            inertia(1., 3.),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn closure(ball: bool) -> BilateralLoop {
    if ball {
        BilateralLoop::PointCoincidence(PointCoincidenceLoop {
            a: MaterialPoint {
                link: 2,
                point_body: Vec3::new(1., 0., 0.),
            },
            b: point(5),
        })
    } else {
        BilateralLoop::Distance(DistanceLoop {
            a: point(1),
            b: point(3),
            length_m: 1.,
        })
    }
}
fn initial(ball: bool, locked: bool) -> PrescribedBaseMotionState {
    let mut q = vec![Q0[0], Q0[1], Q0[0]];
    let mut v = vec![V0[0], if locked { 0. } else { V0[1] }, V0[0]];
    if ball {
        q.extend([0., 0.]);
        v.extend([0., 0.]);
    }
    PrescribedBaseMotionState {
        time_s: T0,
        base: support(T0),
        joint_position: q,
        joint_velocity: v,
    }
}
fn loads(ball: bool) -> SpatialMotionLoads<'static> {
    let n = if ball { 5 } else { 3 };
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &HELD[..n],
            external_body_wrench: &WRENCHES[..n + 1],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    }
}
fn lock(model: &ArticulatedModel, index: usize, position: f64) -> JointConstraint {
    let mut coefficients = vec![0.; model.dof_count()];
    coefficients[index] = 1.;
    JointConstraint {
        coordinate: if index == 1 {
            ConstraintCoordinate::Angle
        } else {
            ConstraintCoordinate::Length
        },
        coefficients,
        position,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }
}
fn diagonal() -> SpatialBushing {
    SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: pose(0.),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: pose(0.),
        },
        rest_translation_a: Vec3::new(1., 0., 0.),
        translation_stiffness: Vec3::new(K[0][0], 0., 0.),
        translation_damping: Vec3::new(C[0][0], 0., 0.),
        rotation_stiffness: Vec3::new(K[1][1], 0., 0.),
        rotation_damping: Vec3::new(C[1][1], 0., 0.),
    }
}
fn coordinate_scales() -> BushingCoordinateScales {
    BushingCoordinateScales {
        translation_m: 0.25,
        rotation_rad: 1.,
    }
}
fn coupled() -> CoupledSpatialBushing {
    let embed = |a: [[f64; 2]; 2]| {
        let mut m = [[0.; 6]; 6];
        for i in 0..2 {
            for j in 0..2 {
                m[3 * i][3 * j] = a[i][j];
            }
        }
        PassiveBushingMatrix::try_new(m, coordinate_scales()).unwrap()
    };
    let b = diagonal();
    CoupledSpatialBushing::try_new(b.a, b.b, b.rest_translation_a, embed(K), embed(C)).unwrap()
}

fn environment(
    model: &ArticulatedModel,
    stage: &PrescribedBaseMotionState,
    wrenches: &mut [Wrench],
    effort: &mut [f64],
) -> Result<f64, PrescribedBaseMotionError> {
    assert!(wrenches.iter().all(|w| *w == Wrench::default()));
    assert_eq!(effort, &HELD[..effort.len()]);
    for i in 0..2 {
        effort[i] -= 2. * stage.joint_position[i] + 0.3 * stage.joint_velocity[i];
    }
    let fk = forward_kinematics(
        model,
        stage.base,
        &stage.joint_position,
        &stage.joint_velocity,
    )?;
    let pose = fk.world_from_link[2];
    let twist = fk.body_twist[2];
    let mut power = 0.;
    // Real off-axis material points: translation and spin both enter drag.
    for sign in [-1., 1.] {
        let arm = Vec3::new(0., sign * 0.3, 0.);
        let point = pose.transform_point(arm)?;
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(arm))?;
        let force = Vec3::new(
            -0.5 * (DRAG[0] * velocity.x + 3. * (point.x - 1.2)),
            -DRAG[1] / 0.18 * velocity.y,
            -DRAG[1] / 0.18 * velocity.z,
        );
        let local = pose.rotation().inverse().rotate(force)?;
        wrenches[2].force = wrenches[2].force + local;
        wrenches[2].torque = wrenches[2].torque + arm.cross(local);
        power += force.dot(velocity);
    }
    Ok(power)
}

fn reference(ball: bool, locked: bool, start: &PrescribedBaseMotionState, h: f64) -> [f64; 2] {
    // Independent reduction: prescribed carrier; one combined translating mass
    // and one rotor inertia. No native mass, mobility or nonlinear solve is used.
    let mass = [if ball { 5. } else { 3. }, 2.];
    let base = path(start.time_s, h).midpoint;
    let speed = [base.twist_body.linear.x, base.twist_body.angular.x];
    let half = 0.5 * h;
    let a: [[f64; 2]; 2] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let diag = if i == j { 1. } else { 0. };
            mass[i] * diag
                + half * (C[i][j] + diag * (DRAG[i] + 0.3))
                + half * half * (K[i][j] + diag * (2. + if i == 0 { 3. } else { 0. }))
        })
    });
    let rhs: [f64; 2] = core::array::from_fn(|i| {
        mass[i] * start.joint_velocity[i]
            + half
                * (HELD[i]
                    - 2. * Q0[i]
                    - K[i][0] * Q0[0]
                    - K[i][1] * Q0[1]
                    - DRAG[i] * speed[i]
                    - mass[i] * [A, ALPHA][i]
                    - if i == 0 {
                        3. * (base.world_from_base.translation().x + 1. + Q0[0] - 1.2)
                    } else {
                        0.
                    })
    });
    if locked {
        [rhs[0] / a[0][0], 0.]
    } else {
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        [
            (a[1][1] * rhs[0] - a[0][1] * rhs[1]) / det,
            (a[0][0] * rhs[1] - a[1][0] * rhs[0]) / det,
        ]
    }
}

fn accounts(s: &PrescribedBilateralLoopStep) {
    let m = &s.midpoint.motion;
    let actual = Wrench::new(
        s.momentum_change_world.torque
            - m.external_angular_impulse_world
            - s.support_impulse_world.torque
            - s.constraint_impulse_world.torque,
        s.momentum_change_world.force
            - m.external_linear_impulse_world
            - s.support_impulse_world.force
            - s.constraint_impulse_world.force,
    );
    vector(actual.force, s.momentum_balance_residual_world.force, 2e-10);
    vector(
        actual.torque,
        s.momentum_balance_residual_world.torque,
        2e-10,
    );
    near(
        s.kinetic_energy_change + s.midpoint.stored_energy_after - s.midpoint.stored_energy_before
            + s.midpoint.dissipated_energy
            - m.actuator_work
            - m.external_work
            - s.support_work
            - s.constraint_impulse_work,
        s.mechanical_work_residual,
        2e-10,
    );
    assert!(s.position_residual_after_m.iter().all(|x| x.abs() < 1e-9));
    assert!(s.velocity_residual_after_mps.iter().all(|x| x.abs() < 1e-9));
    for p in &s.point_residual_after {
        assert!(p.relative_velocity_world_mps.norm() < 1e-9);
    }
}

#[test]
fn closed_dragged_payloads_and_full_bushings_keep_actual_accelerating_spinning_support_exchange() {
    for ball in [false, true] {
        for locked in [false, true] {
            let model = model(ball);
            let start = initial(ball, locked);
            let saved = start.clone();
            let rows = if locked {
                vec![lock(&model, 1, Q0[1])]
            } else {
                Vec::new()
            };
            let h = 0.4;
            let half = h / 2.;
            let expected = reference(ball, locked, &start, h);
            let mut last = None;
            let mut calls = 0;
            let step = step_with_coupled_bushings_and_forcing(
                &model,
                &start,
                path(T0, h),
                loads(ball),
                &[closure(ball)],
                &rows,
                &[coupled()],
                h,
                config(),
                |elapsed, stage, wrenches, effort| {
                    near(elapsed, half, 1e-15);
                    near(stage.time_s, T0 + half, 1e-15);
                    assert_eq!(stage.base, path(T0, h).midpoint);
                    let power = environment(&model, stage, wrenches, effort)?;
                    calls += 1;
                    last = Some((stage.clone(), power, wrenches[2], effort.to_vec()));
                    Ok(())
                },
                || false,
            )
            .unwrap();
            let (accepted, power, wrench, effort) = last.unwrap();
            assert_eq!(accepted, step.midpoint_state);
            assert_eq!(calls, step.residual_evaluations);
            assert_eq!(step.state_after.base, path(T0, h).endpoint);
            near(step.state_after.time_s, T0 + h, 1e-15);
            for i in 0..2 {
                near(accepted.joint_velocity[i], expected[i], 2e-10);
                near(
                    step.state_after.joint_position[i],
                    Q0[i] + h * expected[i],
                    2e-10,
                );
                near(
                    step.state_after.joint_velocity[i],
                    2. * expected[i] - start.joint_velocity[i],
                    2e-10,
                );
            }
            let [vx, spin] = expected;
            let base = path(T0, h).midpoint;
            let u = base.twist_body.linear.x;
            let omega = base.twist_body.angular.x;
            let force = -DRAG[0] * (u + vx)
                - 3. * (base.world_from_base.translation().x + 1. + Q0[0] + half * vx - 1.2);
            let torque = -DRAG[1] * (omega + spin);
            near(wrench.force.x, force, 2e-9);
            near(wrench.torque.x, torque, 2e-9);
            near(power, force * (u + vx) + torque * (omega + spin), 2e-9);
            near(step.midpoint.motion.external_work, h * power, 2e-9);
            near(
                step.midpoint.motion.actuator_work,
                h * (effort[0] * vx + effort[1] * spin),
                2e-9,
            );
            let moving = if ball { 5. } else { 3. };
            let total_spin = if ball { 16. } else { 10. };
            let acceleration = (vx - start.joint_velocity[0]) / half;
            let spin_acceleration = (spin - start.joint_velocity[1]) / half;
            let support_force = (3. + moving) * A + moving * acceleration - force;
            let support_torque = total_spin * ALPHA + 2. * spin_acceleration - torque;
            vector(
                step.support_impulse_world.force,
                Vec3::new(h * support_force, 0., 0.),
                2e-9,
            );
            vector(
                step.support_impulse_world.torque,
                Vec3::new(
                    h * support_torque,
                    0.4 * h * support_force,
                    0.7 * h * support_force,
                ),
                2e-9,
            );
            near(
                step.support_work,
                h * (support_force * u + support_torque * omega),
                2e-9,
            );
            let reaction = half * (if ball { 3. } else { 1. }) * (A + acceleration);
            assert!(reaction.abs() > 1e-3);
            near(step.initial_loop_impulses_ns[0], reaction, 2e-9);
            near(step.final_loop_impulses_ns[0], reaction, 2e-9);
            if locked {
                let q = &accepted.joint_position;
                let bushing_torque = K[1][0] * q[0] + K[1][1] * q[1] + C[1][0] * vx;
                let lambda = half * (2. * ALPHA - effort[1] + bushing_torque - torque);
                near(step.initial_joint_multipliers[0], lambda, 2e-9);
                near(step.final_joint_multipliers[0], lambda, 2e-9);
                near(
                    step.initial_reaction.support_work,
                    lambda * start.base.twist_body.angular.x,
                    2e-9,
                );
                near(
                    step.final_reaction.support_work,
                    lambda * step.state_after.base.twist_body.angular.x,
                    2e-9,
                );
                assert!(step.initial_reaction.support_work.abs() > 1e-3);
                assert!(step.final_reaction.support_work.abs() > 1e-3);
            }
            let kinetic = |state: &PrescribedBaseMotionState| {
                let u = state.base.twist_body.linear.x;
                let w = state.base.twist_body.angular.x;
                let v = state.joint_velocity[0];
                let s = state.joint_velocity[1];
                0.5 * (3. + moving) * u * u
                    + moving * u * v
                    + 0.5 * moving * v * v
                    + 0.5 * total_spin * w * w
                    + 2. * w * s
                    + s * s
            };
            near(
                step.kinetic_energy_change,
                kinetic(&step.state_after) - kinetic(&start),
                2e-9,
            );
            near(step.constraint_impulse_work, 0., 2e-9);
            near(step.mechanical_work_residual, 0., 2e-9);
            vector(
                step.momentum_balance_residual_world.force,
                Vec3::default(),
                2e-9,
            );
            vector(
                step.momentum_balance_residual_world.torque,
                Vec3::default(),
                2e-9,
            );
            accounts(&step);
            assert_eq!(start, saved);
        }
    }
}

#[test]
fn mobile_scalar_root_and_material_lock_keep_native_support_work_and_diagonal_parity() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "mobile root",
            None,
            pose(0.),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(3., 3.),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(1.),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2., 2.),
        ),
    ])
    .unwrap();
    let start = PrescribedBaseMotionState {
        time_s: T0,
        base: support(T0),
        joint_position: vec![0.2, 0.1],
        joint_velocity: vec![0., 0.],
    };
    let loops = [BilateralLoop::Distance(DistanceLoop {
        a: point(0),
        b: point(1),
        length_m: 1.1,
    })];
    let rows = [lock(&model, 0, 0.2)];
    let body = [
        Wrench::default(),
        Wrench::new(Vec3::default(), Vec3::new(2.5, 0., 0.)),
    ];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[7., 1.6],
            external_body_wrench: &body,
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let b = SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: pose(0.),
        },
        b: MaterialFrame {
            link: 1,
            link_from_attachment: pose(0.),
        },
        rest_translation_a: Vec3::new(1., 0., 0.),
        translation_stiffness: Vec3::new(4., 0., 0.),
        translation_damping: Vec3::new(0.4, 0., 0.),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
    };
    let full = CoupledSpatialBushing::try_from_diagonal(b, coordinate_scales()).unwrap();
    let h = 0.3;
    let plain = step_with_bushings(
        &model,
        &start,
        path(T0, h),
        loads,
        &loops,
        &rows,
        &[b],
        h,
        config(),
        || false,
    )
    .unwrap();
    let forced = step_with_bushings_and_forcing(
        &model,
        &start,
        path(T0, h),
        loads,
        &loops,
        &rows,
        &[b],
        h,
        config(),
        |_, _, _, _| Ok(()),
        || false,
    )
    .unwrap();
    let coupled = step_with_coupled_bushings_and_forcing(
        &model,
        &start,
        path(T0, h),
        loads,
        &loops,
        &rows,
        &[full],
        h,
        config(),
        |_, _, _, _| Ok(()),
        || false,
    )
    .unwrap();
    assert_eq!(plain, forced);
    assert_eq!(plain, coupled);
    let mid = path(T0, h).midpoint;
    let force = 5. * A - 2.5;
    near(plain.support_impulse_world.force.x, h * force, 2e-10);
    near(plain.support_impulse_world.torque.x, h * 5. * ALPHA, 2e-10);
    near(
        plain.support_work,
        h * (force * mid.twist_body.linear.x + 5. * ALPHA * mid.twist_body.angular.x),
        2e-10,
    );
    assert!(plain.initial_reaction.support_work.abs() > 1e-4);
    assert!(plain.final_reaction.support_work.abs() > 1e-4);
    for i in 0..2 {
        near(
            plain.state_after.joint_position[i],
            start.joint_position[i],
            2e-11,
        );
        near(plain.state_after.joint_velocity[i], 0., 2e-11);
    }
    near(plain.constraint_impulse_work, 0., 2e-10);
    near(plain.mechanical_work_residual, 0., 2e-10);
    accounts(&plain);
}

#[test]
fn supported_loop_trial_final_sample_and_budget_refuse_without_publishing_state() {
    let model = model(true);
    let start = initial(true, false);
    let saved = start.clone();
    let loops = [closure(true)];
    let bushings = [coupled()];
    let h = 0.4;
    let run = || {
        step_with_coupled_bushings_and_forcing(
            &model,
            &start,
            path(T0, h),
            loads(true),
            &loops,
            &[],
            &bushings,
            h,
            config(),
            |_, s, w, e| environment(&model, s, w, e).map(|_| ()),
            || false,
        )
        .unwrap()
    };
    let reference = run();
    let last = reference.residual_evaluations;
    for mode in 0..6 {
        let calls = Cell::new(0);
        let stop = Cell::new(false);
        let fail = if mode == 0 { 2 } else { last };
        let controls = if mode == 5 {
            BilateralLoopConfig {
                nonlinear: ImplicitSpatialConfig {
                    maximum_evaluations: last - 1,
                    ..config().nonlinear
                },
                ..config()
            }
        } else {
            config()
        };
        let result = step_with_coupled_bushings_and_forcing(
            &model,
            &start,
            path(T0, h),
            loads(true),
            &loops,
            &[],
            &bushings,
            h,
            controls,
            |_, s, w, e| {
                calls.set(calls.get() + 1);
                environment(&model, s, w, e)?;
                if calls.get() == fail {
                    match mode {
                        0 | 1 => {
                            return Err(PrescribedBaseMotionError::Invalid("environment refused"));
                        }
                        2 => w[2].force.x = f64::NAN,
                        3 => e[0] = f64::INFINITY,
                        4 => stop.set(true),
                        _ => {}
                    }
                }
                Ok(())
            },
            || stop.get(),
        );
        assert!(result.is_err(), "mode {mode}");
        if mode <= 1 {
            assert_eq!(
                result,
                Err(BilateralLoopError::SupportedMotion(
                    ImplicitPrescribedSpatialError::Motion(PrescribedBaseMotionError::Invalid(
                        "environment refused"
                    ))
                ))
            );
            assert_eq!(calls.get(), fail);
        }
        if mode == 5 {
            assert_eq!(
                result,
                Err(BilateralLoopError::Motion(
                    ImplicitSpatialError::EvaluationBudget {
                        used: last - 1,
                        maximum: last - 1
                    }
                ))
            );
        }
        assert_eq!(start, saved);
    }
    let bad = [BilateralLoop::Distance(DistanceLoop {
        a: point(1),
        b: point(3),
        length_m: 1.01,
    })];
    assert!(matches!(
        step_with_bushings(
            &model,
            &start,
            path(T0, h),
            loads(true),
            &bad,
            &[],
            &[],
            h,
            config(),
            || false
        ),
        Err(BilateralLoopError::Position {
            stage: "initial",
            ..
        })
    ));
    // Permissive numerical chart probes must not waive the original physical
    // endpoint update admission in the shared supported midpoint owner.
    let mobile = ArticulatedModel::new(vec![Link::new(
        "unrepresentable physical travel",
        None,
        pose(0.),
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
        inertia(1., 1.),
    )])
    .unwrap();
    let tiny = PrescribedBaseMotionState {
        time_s: T0,
        base: BaseState::stationary(Se3::identity()),
        joint_position: vec![1.],
        joint_velocity: vec![1e-200],
    };
    let none = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.],
            external_body_wrench: &WRENCHES[..1],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    assert_eq!(
        step_prescribed_base_bushings_implicit_with_forcing(
            &mobile,
            &tiny,
            PrescribedBasePath::stationary(Se3::identity()),
            none,
            &[],
            0.1,
            config().nonlinear,
            |_, _, _, _| Ok(()),
            || false,
        ),
        Err(ImplicitPrescribedSpatialError::Motion(
            PrescribedBaseMotionError::Invalid("unrepresentable joint update")
        ))
    );
    assert_eq!(run(), reference);
}

fn rotating_support(time: f64) -> BaseState {
    let t = time - T0;
    let axis = Vec3::new(0.2, -0.3, 0.4).scale(1. / 0.29_f64.sqrt());
    let r = rotation(Vec3::new(0.3, -0.2, 0.1))
        .compose(rotation(axis.scale(0.7 * t + 0.1 * t * t)))
        .unwrap();
    let a = Vec3::new(0.15, -0.09, 0.07);
    let u0 = Vec3::new(0.2, 0.1, -0.1);
    let x = Vec3::new(0.2, -0.3, 0.4) + u0.scale(t) + a.scale(0.5 * t * t);
    let v = r.inverse().rotate(u0 + a.scale(t)).unwrap();
    let w = axis.scale(0.7 + 0.2 * t);
    BaseState::prescribed(
        Se3::from_parts(r, x).unwrap(),
        Twist::new(w, v),
        Twist::new(axis.scale(0.2), r.inverse().rotate(a).unwrap() - w.cross(v)),
    )
}
fn rotating_path(t: f64, h: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: rotating_support(t + 0.5 * h),
        endpoint: rotating_support(t + h),
    }
}
fn triangle(reduced: bool) -> ArticulatedModel {
    let root = Link::new(
        "driven carrier",
        None,
        pose(0.),
        JointModel::FIXED,
        SpatialInertia::new(
            3.,
            Vec3::default(),
            Mat3 {
                m: [2., 0.2, 0.1, 0.2, 2.5, 0.15, 0.1, 0.15, 3.],
            },
        )
        .unwrap(),
    );
    let hinge = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
    if reduced {
        // Independently combine the original arms at +pi/2, retaining the full
        // parallel-axis tensor and off-diagonal entries in arm A's body frame.
        let compound = SpatialInertia::new(
            2.,
            Vec3::new(0.25, 0.25, 0.),
            Mat3 {
                m: [0.425, 0.125, 0., 0.125, 0.425, 0., 0., 0., 0.75],
            },
        )
        .unwrap();
        ArticulatedModel::new(vec![
            root,
            Link::new("compound", Some(0), pose(0.), hinge, compound),
        ])
        .unwrap()
    } else {
        let rod = SpatialInertia::new(
            1.,
            Vec3::new(0.5, 0., 0.),
            Mat3 {
                m: [0.1, 0., 0., 0., 0.2, 0., 0., 0., 0.25],
            },
        )
        .unwrap();
        ArticulatedModel::new(vec![
            root,
            Link::new("arm A", Some(0), pose(0.), hinge, rod),
            Link::new("arm B", Some(0), pose(0.), hinge, rod),
        ])
        .unwrap()
    }
}
fn triangle_environment(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    reduced: bool,
    wrenches: &mut [Wrench],
    effort: &mut [f64],
) -> Result<(), PrescribedBaseMotionError> {
    let fk = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    let points = [
        (1, Vec3::new(0.7, -0.1, 0.2)),
        if reduced {
            // Arm B's actual material point transformed into the compound frame.
            (1, Vec3::new(-0.2, 0.9, -0.1))
        } else {
            (2, Vec3::new(0.9, 0.2, -0.1))
        },
    ];
    let time = state.time_s - T0;
    for (link, arm) in points {
        let p = fk.world_from_link[link];
        let twist = fk.body_twist[link];
        let point = p.transform_point(arm)?;
        let velocity = p
            .rotation()
            .rotate(twist.linear + twist.angular.cross(arm))?;
        let relative = Vec3::new(0.08, -0.03, 0.02) - velocity;
        let force = relative.scale(0.5 + 0.8 * relative.norm())
            + (Vec3::new(0.5, -0.2, 0.1) - point).scale(0.2)
            + Vec3::new(0.15 + 0.1 * time, -0.07, 0.03 - 0.05 * time);
        let local = p.rotation().inverse().rotate(force)?;
        wrenches[link].force = wrenches[link].force + local;
        wrenches[link].torque = wrenches[link].torque + arm.cross(local);
    }
    if reduced {
        effort[0] -= 0.4 * state.joint_velocity[0];
    } else {
        effort[0] -= 0.3 * state.joint_velocity[0];
        effort[1] -= 0.1 * state.joint_velocity[1];
    }
    Ok(())
}

#[test]
fn nonlinear_closed_triangle_on_rotating_accelerating_support_converges_to_full_inertia_reduction()
{
    let model = triangle(false);
    let reduced = triangle(true);
    let half_pi = core::f64::consts::FRAC_PI_2;
    let initial = PrescribedBaseMotionState {
        time_s: T0,
        base: rotating_support(T0),
        joint_position: vec![0.4, 0.4 + half_pi],
        joint_velocity: vec![0.7, 0.7],
    };
    let mut reference = PrescribedBaseMotionState {
        joint_position: vec![0.4],
        joint_velocity: vec![0.7],
        ..initial.clone()
    };
    let closure = BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::new(1., 0., 0.),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::new(1., 0., 0.),
        },
        length_m: 2.0_f64.sqrt(),
    });
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.8, -0.2],
            external_body_wrench: &WRENCHES[..3],
            gravity_world: Vec3::new(0., -0.4, 0.),
        },
        world_points: &[],
        connections: &[],
    };
    let reduced_loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.6],
            external_body_wrench: &WRENCHES[..2],
            ..loads.held
        },
        ..loads
    };
    let duration = 0.3;
    let mut reference_support = Vec3::default();
    let mut reference_couple = Vec3::default();
    let mut reference_work = 0.;
    for _ in 0..256 {
        let h = duration / 256.;
        let result = step_prescribed_base_bushings_implicit_with_forcing(
            &reduced,
            &reference,
            rotating_path(reference.time_s, h),
            reduced_loads,
            &[],
            h,
            config().nonlinear,
            |elapsed, s, w, e| {
                near(s.time_s, reference.time_s + elapsed, 2e-15);
                triangle_environment(&reduced, s, true, w, e)
            },
            || false,
        )
        .unwrap();
        reference_support = reference_support + result.spatial.motion.support_linear_impulse_world;
        reference_couple = reference_couple + result.spatial.motion.support_angular_impulse_world;
        reference_work += result.spatial.motion.support_work;
        reference = result.spatial.motion.state_after;
    }
    let mut errors = Vec::new();
    let mut support_errors = Vec::new();
    let mut work_errors = Vec::new();
    let mut momentum_defects = Vec::new();
    for count in [8, 16, 32] {
        let mut state = initial.clone();
        let mut support = Vec3::default();
        let mut couple = Vec3::default();
        let mut work = 0.;
        let mut reactions = 0.;
        let mut defect = 0.;
        for _ in 0..count {
            let h = duration / f64::from(count);
            let mut calls = 0;
            let mut accepted = None;
            let result = step_with_bushings_and_forcing(
                &model,
                &state,
                rotating_path(state.time_s, h),
                loads,
                &[closure],
                &[],
                &[],
                h,
                config(),
                |elapsed, s, w, e| {
                    near(s.time_s, state.time_s + elapsed, 2e-15);
                    calls += 1;
                    accepted = Some(s.clone());
                    triangle_environment(&model, s, false, w, e)
                },
                || false,
            )
            .unwrap();
            assert_eq!(calls, result.residual_evaluations);
            assert_eq!(accepted, Some(result.midpoint_state.clone()));
            accounts(&result);
            reactions +=
                result.initial_loop_impulses_ns[0].abs() + result.final_loop_impulses_ns[0].abs();
            support = support + result.support_impulse_world.force;
            couple = couple + result.support_impulse_world.torque;
            work += result.support_work;
            defect += result.momentum_balance_residual_world.force.norm()
                + result.momentum_balance_residual_world.torque.norm();
            state = result.state_after;
            near(
                state.joint_position[1] - state.joint_position[0],
                half_pi,
                1e-9,
            );
            near(state.joint_velocity[1], state.joint_velocity[0], 1e-9);
            let fk = forward_kinematics(
                &model,
                state.base,
                &state.joint_position,
                &state.joint_velocity,
            )
            .unwrap();
            let a = fk.world_from_link[1]
                .transform_point(Vec3::new(1., 0., 0.))
                .unwrap();
            let b = fk.world_from_link[2]
                .transform_point(Vec3::new(1., 0., 0.))
                .unwrap();
            near((b - a).norm(), 2.0_f64.sqrt(), 1e-9);
        }
        assert!(reactions > 1e-3);
        errors.push(
            ((state.joint_position[0] - reference.joint_position[0]).powi(2)
                + (state.joint_velocity[0] - reference.joint_velocity[0]).powi(2))
            .sqrt(),
        );
        support_errors
            .push((support - reference_support).norm() + (couple - reference_couple).norm());
        work_errors.push((work - reference_work).abs());
        momentum_defects.push(defect);
    }
    eprintln!(
        "supported triangle state errors {errors:?}; support impulse errors {support_errors:?}; support work errors {work_errors:?}; native momentum defects {momentum_defects:?}"
    );
    for sequence in [&errors, &support_errors, &work_errors] {
        for pair in sequence.windows(2) {
            assert!((3.5..4.5).contains(&(pair[0] / pair[1])), "{sequence:?}");
        }
    }
    assert!(errors[2] < 1e-4);
    // The native prescribed-path midpoint has a finite integration defect. It is
    // retained and decreases with refinement, rather than defining the support.
    assert!(momentum_defects[0] > 1e-8);
    assert!(momentum_defects[2] < momentum_defects[0] / 12.);
}
