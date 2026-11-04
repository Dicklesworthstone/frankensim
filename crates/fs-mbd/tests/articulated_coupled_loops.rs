//! Full passive translation/rotation coupling inside native closed mechanisms.
//!
//! G1 reduced two-mode midpoint oracles retain both carrier recoils, rigid
//! reactions, finite storage and loss. G0/G4 retain parity and atomic refusal.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, BilateralLoopStep, DistanceLoop,
    PointCoincidenceLoop, joint, step_free_floating_bilateral_loops,
    step_free_floating_bilateral_loops_with_bushings,
    step_free_floating_bilateral_loops_with_coupled_bushings,
    step_free_floating_bilateral_loops_with_coupled_bushings_and_effort,
    step_free_floating_bilateral_loops_with_coupled_bushings_and_forcing,
    step_free_floating_bilateral_loops_with_effort,
    step_free_floating_bilateral_loops_with_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics,
};

const ZERO: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};
const WRENCHES: [Wrench; 6] = [Wrench::new(ZERO, ZERO); 6];
const HELD: [f64; 5] = [0.3, -0.2, 0.0, 0.0, 0.0];
const Q0: [f64; 2] = [0.3, -0.2];
const V0: [f64; 2] = [-0.1, 0.4];
const K: [[f64; 2]; 2] = [[12.0, 4.0], [4.0, 8.0]];
const C: [[f64; 2]; 2] = [[0.6, 0.1], [0.1, 0.4]];

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8, "{a:.16e} != {b:.16e}");
}

fn pose(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, 0.0)).unwrap()
}

fn model(point: bool) -> ArticulatedModel {
    let inertia = |mass, spin| {
        SpatialInertia::new(
            mass,
            ZERO,
            Mat3 {
                m: [spin, 0.0, 0.0, 0.0, spin, 0.0, 0.0, 0.0, spin],
            },
        )
        .unwrap()
    };
    let slide = |axis| JointModel::prismatic(axis, None).unwrap();
    let mut links = vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, 3.0),
        ),
        Link::new(
            "slider",
            Some(0),
            pose(1.0),
            slide(Vec3::new(1.0, 0.0, 0.0)),
            inertia(1.0, 2.0),
        ),
        Link::new(
            "rotor",
            Some(1),
            Se3::identity(),
            JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
        Link::new(
            "follower x",
            Some(0),
            pose(2.0),
            slide(Vec3::new(1.0, 0.0, 0.0)),
            inertia(1.0, 3.0),
        ),
    ];
    if point {
        // The two extra physical sliders give the ball closure three independent
        // reaction directions. Their masses/inertias remain in the native tree.
        links.push(Link::new(
            "follower y",
            Some(3),
            Se3::identity(),
            slide(Vec3::new(0.0, 1.0, 0.0)),
            inertia(1.0, 3.0),
        ));
        links.push(Link::new(
            "follower z",
            Some(4),
            Se3::identity(),
            slide(Vec3::new(0.0, 0.0, 1.0)),
            inertia(1.0, 3.0),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}

fn closure(point: bool) -> BilateralLoop {
    if point {
        BilateralLoop::PointCoincidence(PointCoincidenceLoop {
            a: MaterialPoint {
                link: 2,
                point_body: Vec3::new(1.0, 0.0, 0.0),
            },
            b: MaterialPoint {
                link: 5,
                point_body: ZERO,
            },
        })
    } else {
        BilateralLoop::Distance(DistanceLoop {
            a: MaterialPoint {
                link: 1,
                point_body: ZERO,
            },
            b: MaterialPoint {
                link: 3,
                point_body: ZERO,
            },
            length_m: 1.0,
        })
    }
}

fn mass(point: bool) -> [f64; 2] {
    // Independent reduction: three/five translating kg against the 3 kg
    // carrier; rotor spin 2 against the other bodies' spin 8/14.
    if point {
        [15.0 / 8.0, 7.0 / 4.0]
    } else {
        [3.0 / 2.0, 8.0 / 5.0]
    }
}

fn recoil(point: bool) -> [f64; 2] {
    if point {
        [5.0 / 8.0, 1.0 / 8.0]
    } else {
        [0.5, 0.2]
    }
}

fn initial(point: bool) -> FreeFloatingMotionState {
    let fraction = recoil(point);
    let mut q = vec![Q0[0], Q0[1], Q0[0]];
    let mut v = vec![V0[0], V0[1], V0[0]];
    if point {
        q.extend([0.0, 0.0]);
        v.extend([0.0, 0.0]);
    }
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(
                Vec3::new(-fraction[1] * V0[1], 0.0, 0.0),
                Vec3::new(-fraction[0] * V0[0], 0.0, 0.0),
            ),
        ),
        joint_position: q,
        joint_velocity: v,
    }
}

fn loads(point: bool) -> SpatialMotionLoads<'static> {
    let dofs = if point { 5 } else { 3 };
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &HELD[..dofs],
            external_body_wrench: &WRENCHES[..dofs + 1],
            gravity_world: ZERO,
        },
        connections: &[],
        world_points: &[],
    }
}

fn diagonal() -> SpatialBushing {
    SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::new(1.0, 0.0, 0.0),
        translation_stiffness: Vec3::new(K[0][0], 0.0, 0.0),
        translation_damping: Vec3::new(C[0][0], 0.0, 0.0),
        rotation_stiffness: Vec3::new(K[1][1], 0.0, 0.0),
        rotation_damping: Vec3::new(C[1][1], 0.0, 0.0),
    }
}

fn scales() -> BushingCoordinateScales {
    BushingCoordinateScales {
        translation_m: 0.25,
        rotation_rad: 1.0,
    }
}

fn coupled() -> CoupledSpatialBushing {
    let embed = |a: [[f64; 2]; 2]| {
        let mut matrix = [[0.0; 6]; 6];
        for i in 0..2 {
            for j in 0..2 {
                matrix[3 * i][3 * j] = a[i][j];
            }
        }
        PassiveBushingMatrix::try_new(matrix, scales()).unwrap()
    };
    let b = diagonal();
    CoupledSpatialBushing::try_new(b.a, b.b, b.rest_translation_a, embed(K), embed(C)).unwrap()
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

fn control(
    _: f64,
    stage: &FreeFloatingMotionState,
    output: &mut [f64],
) -> Result<(), FreeFloatingMotionError> {
    assert_eq!(output, &HELD[..output.len()]);
    for i in 0..2 {
        output[i] -= 2.0 * stage.joint_position[i] + 0.3 * stage.joint_velocity[i];
    }
    Ok(())
}

fn midpoint(point: bool, feedback: bool, h: f64) -> [f64; 2] {
    // Independent 2x2 implicit midpoint solve after eliminating the closed
    // mechanism and both conserved momenta; no native mobility is used here.
    let m = mass(point);
    let a: [[f64; 2]; 2] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let delta = if i == j { 1.0 } else { 0.0 };
            let active = if feedback { delta } else { 0.0 };
            delta
                + h * (C[i][j] + 0.3 * active) / (2.0 * m[i])
                + h * h * (K[i][j] + 2.0 * active) / (4.0 * m[i])
        })
    });
    let rhs: [f64; 2] = core::array::from_fn(|i| {
        let kq: f64 = (0..2).map(|j| K[i][j] * Q0[j]).sum();
        V0[i] + h * (HELD[i] - kq - if feedback { 2.0 * Q0[i] } else { 0.0 }) / (2.0 * m[i])
    });
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    [
        (a[1][1] * rhs[0] - a[0][1] * rhs[1]) / det,
        (a[0][0] * rhs[1] - a[1][0] * rhs[0]) / det,
    ]
}

fn check(step: &BilateralLoopStep, point: bool, feedback: bool, h: f64) {
    let vm = midpoint(point, feedback, h);
    let qm: [f64; 2] = core::array::from_fn(|i| Q0[i] + 0.5 * h * vm[i]);
    let q1: [f64; 2] = core::array::from_fn(|i| Q0[i] + h * vm[i]);
    let v1: [f64; 2] = core::array::from_fn(|i| 2.0 * vm[i] - V0[i]);
    let u: [f64; 2] = core::array::from_fn(|i| {
        HELD[i]
            - if feedback {
                2.0 * qm[i] + 0.3 * vm[i]
            } else {
                0.0
            }
    });
    for i in 0..2 {
        close(step.state_after.joint_position[i], q1[i]);
        close(step.state_after.joint_velocity[i], v1[i]);
    }
    close(step.state_after.joint_position[2], q1[0]);
    close(step.state_after.joint_velocity[2], v1[0]);
    for i in 3..step.state_after.joint_position.len() {
        close(step.state_after.joint_position[i], 0.0);
        close(step.state_after.joint_velocity[i], 0.0);
    }
    let fraction = recoil(point);
    close(
        step.state_after.base.twist_body.linear.x,
        -fraction[0] * v1[0],
    );
    close(
        step.state_after.base.twist_body.angular.x,
        -fraction[1] * v1[1],
    );
    close(
        step.state_after.base.world_from_base.translation().x,
        -fraction[0] * h * vm[0],
    );
    let momentum = step.final_reaction.momentum_after;
    close(momentum.linear_momentum_world.norm(), 0.0);
    close(momentum.angular_momentum_world.norm(), 0.0);
    close(
        step.midpoint.motion.center_of_mass_world_after.x,
        if point {
            1.0 + 5.0 * Q0[0] / 8.0
        } else {
            2.0 / 3.0 + Q0[0] / 2.0
        },
    );
    let m = mass(point);
    close(
        momentum.kinetic_energy,
        0.5 * (m[0] * v1[0] * v1[0] + m[1] * v1[1] * v1[1]),
    );
    let quadratic = |a: [[f64; 2]; 2], x: [f64; 2]| -> f64 {
        (0..2)
            .map(|i| (0..2).map(|j| x[i] * a[i][j] * x[j]).sum::<f64>())
            .sum()
    };
    close(step.midpoint.stored_energy_before, 0.5 * quadratic(K, Q0));
    close(step.midpoint.stored_energy_after, 0.5 * quadratic(K, q1));
    close(step.midpoint.dissipated_energy, h * quadratic(C, vm));
    assert!(step.midpoint.stored_energy_after > 0.0);
    assert!(step.midpoint.dissipated_energy > 0.0);
    close(
        step.midpoint.motion.actuator_work,
        h * (u[0] * vm[0] + u[1] * vm[1]),
    );
    let bushing_force: f64 = (0..2).map(|j| K[0][j] * qm[j] + C[0][j] * vm[j]).sum();
    // Follower mass is 1/3 of the distance payload, or 3/5 of the ball
    // payload. Each original constraint event supplies half its step impulse.
    let reaction = h * (u[0] - bushing_force) * if point { 3.0 / 10.0 } else { 1.0 / 6.0 };
    close(step.initial_impulse_ns[0], reaction);
    close(step.final_impulse_ns[0], reaction);
    for j in 1..step.rows.len() {
        close(step.initial_impulse_ns[j], 0.0);
        close(step.final_impulse_ns[j], 0.0);
    }
    for residual in step
        .position_residual_after_m
        .iter()
        .chain(&step.velocity_residual_after_mps)
    {
        close(*residual, 0.0);
    }
    for residual in &step.point_residual_after {
        close(residual.relative_velocity_world_mps.norm(), 0.0);
    }
    for value in [
        step.reaction_linear_impulse_world,
        step.reaction_angular_impulse_world,
        step.linear_momentum_balance_residual_world,
        step.angular_momentum_balance_residual_world,
    ] {
        close(value.norm(), 0.0);
    }
    close(step.constraint_impulse_work, 0.0);
    close(step.mechanical_work_residual, 0.0);
    close(
        step.mechanical_work_residual,
        step.constraint_impulse_work
            + step.midpoint.mechanical_work_residual
            + step.initial_reaction.kinetic_work_residual
            + step.final_reaction.kinetic_work_residual,
    );
}

#[test]
fn full_passive_cross_coupling_drives_distance_and_ball_closures_with_native_recoil() {
    for point in [false, true] {
        let model = model(point);
        let initial = initial(point);
        let saved = initial.clone();
        let h = 0.25;
        let step = step_free_floating_bilateral_loops_with_coupled_bushings(
            &model,
            &initial,
            loads(point),
            &[closure(point)],
            &[coupled()],
            h,
            config(),
            || false,
        )
        .unwrap();
        check(&step, point, false, h);
        assert_eq!(initial, saved);
    }
}

#[test]
fn stiff_full_matrix_feedback_and_loop_reactions_share_the_independent_midpoint() {
    for point in [false, true] {
        let model = model(point);
        let initial = initial(point);
        let h = 3.0;
        let calls = Cell::new(0);
        let accepted = Cell::new([0.0; 4]);
        let step = step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
            &model,
            &initial,
            loads(point),
            &[closure(point)],
            &[coupled()],
            h,
            config(),
            |elapsed, mid, output| {
                close(elapsed, h / 2.0);
                calls.set(calls.get() + 1);
                accepted.set([
                    mid.joint_position[0],
                    mid.joint_position[1],
                    mid.joint_velocity[0],
                    mid.joint_velocity[1],
                ]);
                control(elapsed, mid, output)
            },
            || false,
        )
        .unwrap();
        assert!(calls.get() > 2);
        assert_eq!(calls.get(), step.residual_evaluations);
        let vm = midpoint(point, true, h);
        for i in 0..2 {
            close(accepted.get()[i], Q0[i] + 0.5 * h * vm[i]);
            close(accepted.get()[i + 2], vm[i]);
        }
        check(&step, point, true, h);
    }
}

#[test]
fn full_diagonal_and_empty_bushings_retain_original_loop_results_exactly() {
    for point in [false, true] {
        let model = model(point);
        let initial = initial(point);
        let loops = [closure(point)];
        let b = diagonal();
        let full = CoupledSpatialBushing::try_from_diagonal(b, scales()).unwrap();
        let h = 0.04;
        assert_eq!(
            step_free_floating_bilateral_loops_with_bushings(
                &model,
                &initial,
                loads(point),
                &loops,
                &[b],
                h,
                config(),
                || false,
            )
            .unwrap(),
            step_free_floating_bilateral_loops_with_coupled_bushings(
                &model,
                &initial,
                loads(point),
                &loops,
                std::slice::from_ref(&full),
                h,
                config(),
                || false,
            )
            .unwrap(),
        );
        assert_eq!(
            step_free_floating_bilateral_loops_with_effort(
                &model,
                &initial,
                loads(point),
                &loops,
                &[b],
                h,
                config(),
                control,
                || false,
            )
            .unwrap(),
            step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
                &model,
                &initial,
                loads(point),
                &loops,
                &[full],
                h,
                config(),
                control,
                || false,
            )
            .unwrap(),
        );
        assert_eq!(
            step_free_floating_bilateral_loops(
                &model,
                &initial,
                loads(point),
                &loops,
                h,
                config(),
                || false,
            )
            .unwrap(),
            step_free_floating_bilateral_loops_with_coupled_bushings(
                &model,
                &initial,
                loads(point),
                &loops,
                &[],
                h,
                config(),
                || false,
            )
            .unwrap(),
        );
        assert_eq!(
            step_free_floating_bilateral_loops_with_effort(
                &model,
                &initial,
                loads(point),
                &loops,
                &[],
                h,
                config(),
                control,
                || false,
            )
            .unwrap(),
            step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
                &model,
                &initial,
                loads(point),
                &loops,
                &[],
                h,
                config(),
                control,
                || false,
            )
            .unwrap(),
        );
    }
}

#[test]
fn full_matrix_loop_feedback_refusal_and_late_cancellation_leave_retry_unchanged() {
    let model = model(true);
    let initial = initial(true);
    let saved = initial.clone();
    let loops = [closure(true)];
    let bushings = [coupled()];
    let run = || {
        step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
            &model,
            &initial,
            loads(true),
            &loops,
            &bushings,
            0.25,
            config(),
            control,
            || false,
        )
        .unwrap()
    };
    let reference = run();
    for failure_call in [2, reference.residual_evaluations] {
        let calls = Cell::new(0);
        let refused = step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
            &model,
            &initial,
            loads(true),
            &loops,
            &bushings,
            0.25,
            config(),
            |elapsed, stage, output| {
                calls.set(calls.get() + 1);
                if calls.get() == failure_call {
                    return Err(FreeFloatingMotionError::InvalidLoad(
                        "coupled loop feedback refused",
                    ));
                }
                control(elapsed, stage, output)
            },
            || false,
        );
        assert_eq!(
            refused,
            Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::InvalidLoad("coupled loop feedback refused")
            )))
        );
        assert_eq!(calls.get(), failure_call);
    }
    let calls = Cell::new(0);
    let stop = Cell::new(false);
    let cancelled = step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
        &model,
        &initial,
        loads(true),
        &loops,
        &bushings,
        0.25,
        config(),
        |elapsed, stage, output| {
            calls.set(calls.get() + 1);
            control(elapsed, stage, output)?;
            stop.set(calls.get() == reference.residual_evaluations);
            Ok(())
        },
        || stop.get(),
    );
    assert_eq!(
        cancelled,
        Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Cancelled
        )))
    );
    assert_eq!(calls.get(), reference.residual_evaluations);
    let budget = BilateralLoopConfig {
        nonlinear: ImplicitSpatialConfig {
            maximum_evaluations: 2,
            ..config().nonlinear
        },
        ..config()
    };
    assert!(matches!(
        step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
            &model,
            &initial,
            loads(true),
            &loops,
            &bushings,
            0.25,
            budget,
            control,
            || false,
        ),
        Err(BilateralLoopError::Motion(
            ImplicitSpatialError::EvaluationBudget {
                used: 2,
                maximum: 2
            }
        ))
    ));
    assert_eq!(initial, saved);
    assert_eq!(run(), reference);
}

const AXIAL_DRAG: f64 = 70.0;
const SPIN_DRAG: f64 = 45.0;
const ENVIRONMENT_K: f64 = 3.0;
const ENVIRONMENT_TARGET_X: f64 = 1.2;

fn environment(
    model: &ArticulatedModel,
    elapsed: f64,
    stage: &FreeFloatingMotionState,
    wrenches: &mut [Wrench],
    efforts: &mut [f64],
) -> Result<f64, FreeFloatingMotionError> {
    assert!(wrenches.iter().all(|w| *w == Wrench::default()));
    control(elapsed, stage, efforts)?;
    let motion = forward_kinematics(
        model,
        BaseState::prescribed(
            stage.base.world_from_base,
            stage.base.twist_body,
            Twist::zero(),
        ),
        &stage.joint_position,
        &stage.joint_velocity,
    )?;
    let pose = motion.world_from_link[2];
    let twist = motion.body_twist[2];
    // Two actual material points straddle the shaft. Their transverse drags
    // give a couple; neither the point spin speed nor its arm is discarded.
    let radius = 0.3;
    let transverse = SPIN_DRAG / (2.0 * radius * radius);
    let mut power = 0.0;
    for sign in [-1.0, 1.0] {
        let arm = Vec3::new(0.0, sign * radius, 0.0);
        let point = pose.transform_point(arm)?;
        let local_velocity = twist.linear + twist.angular.cross(arm);
        let velocity = pose.rotation().rotate(local_velocity)?;
        let force = Vec3::new(
            -0.5 * (AXIAL_DRAG * velocity.x + ENVIRONMENT_K * (point.x - ENVIRONMENT_TARGET_X)),
            -transverse * velocity.y,
            -transverse * velocity.z,
        );
        let local_force = pose.rotation().inverse().rotate(force)?;
        wrenches[2].force = wrenches[2].force + local_force;
        wrenches[2].torque = wrenches[2].torque + arm.cross(local_force);
        power += force.dot(velocity);
    }
    Ok(power)
}

fn environment_midpoint(
    point: bool,
    locked: bool,
    state: &FreeFloatingMotionState,
    h: f64,
) -> [f64; 4] {
    // Independently reduced coordinates [carrier x, carrier angle, payload x,
    // rotor angle]. The loop combines the payloads before this tiny solve.
    let moving = if point { 5.0 } else { 3.0 };
    let spin = if point { 16.0 } else { 10.0 };
    let mass = [
        [3.0 + moving, 0.0, moving, 0.0],
        [0.0, spin, 0.0, 2.0],
        [moving, 0.0, moving, 0.0],
        [0.0, 2.0, 0.0, 2.0],
    ];
    let axial = [1.0, 0.0, 1.0, 0.0];
    let rotation = [0.0, 1.0, 0.0, 1.0];
    let old = [
        state.base.twist_body.linear.x,
        state.base.twist_body.angular.x,
        state.joint_velocity[0],
        state.joint_velocity[1],
    ];
    let half = h / 2.0;
    let mut a: [[f64; 4]; 4] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut stiffness = ENVIRONMENT_K * axial[i] * axial[j];
            let mut damping =
                AXIAL_DRAG * axial[i] * axial[j] + SPIN_DRAG * rotation[i] * rotation[j];
            if i >= 2 && j >= 2 {
                stiffness += K[i - 2][j - 2] + if i == j { 2.0 } else { 0.0 };
                damping += C[i - 2][j - 2] + if i == j { 0.3 } else { 0.0 };
            }
            mass[i][j] + half * damping + half * half * stiffness
        })
    });
    let gap = state.base.world_from_base.translation().x + 1.0 + state.joint_position[0]
        - ENVIRONMENT_TARGET_X;
    let mut rhs: [f64; 4] = core::array::from_fn(|i| {
        let mut source = -ENVIRONMENT_K * gap * axial[i];
        if i >= 2 {
            source += HELD[i - 2]
                - 2.0 * state.joint_position[i - 2]
                - (0..2)
                    .map(|j| K[i - 2][j] * state.joint_position[j])
                    .sum::<f64>();
        }
        (0..4).map(|j| mass[i][j] * old[j]).sum::<f64>() + half * source
    });
    if locked {
        // This independently eliminates the known stationary rotor coordinate.
        a[3] = [0.0, 0.0, 0.0, 1.0];
        rhs[3] = 0.0;
    }
    // Independent small Gaussian solve; no native mass, mobility or nonlinear
    // solver action supplies this oracle's coefficients or arithmetic.
    for k in 0..4 {
        let pivot = (k..4)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, pivot);
        rhs.swap(k, pivot);
        for i in k + 1..4 {
            let multiplier = a[i][k] / a[k][k];
            for j in k + 1..4 {
                a[i][j] -= multiplier * a[k][j];
            }
            rhs[i] -= multiplier * rhs[k];
        }
    }
    let mut result = [0.0; 4];
    for i in (0..4).rev() {
        result[i] = (rhs[i] - (i + 1..4).map(|j| a[i][j] * result[j]).sum::<f64>()) / a[i][i];
    }
    result
}

#[test]
fn material_point_drag_and_full_bushings_drive_closed_payloads_with_both_carrier_recoils() {
    for point in [false, true] {
        let model = model(point);
        for locked in [false, true] {
            let mut start = initial(point);
            start.base.world_from_base = Se3::from_parts(
                So3::exp(So3Tangent::new(Vec3::new(0.2, 0.0, 0.0))).unwrap(),
                Vec3::new(2.0, -0.7, 0.4),
            )
            .unwrap();
            start.base.twist_body.linear.x += 0.4;
            start.base.twist_body.angular.x += 0.3;
            let saved = start.clone();
            let rows = if locked {
                start.joint_velocity[1] = 0.0;
                let mut coefficients = vec![0.0; model.dof_count()];
                coefficients[1] = 1.0;
                vec![JointConstraint {
                    coordinate: ConstraintCoordinate::Angle,
                    coefficients,
                    position: start.joint_position[1],
                    position_tolerance: 1e-9,
                    velocity_tolerance: 1e-9,
                }]
            } else {
                Vec::new()
            };
            let h = 0.4;
            let half = h / 2.0;
            let expected = environment_midpoint(point, locked, &start, h);
            let [u, omega, vx, spin] = expected;
            let qm = [Q0[0] + half * vx, Q0[1] + half * spin];
            let force = -AXIAL_DRAG * (u + vx)
                - ENVIRONMENT_K * (3.0 + qm[0] + half * u - ENVIRONMENT_TARGET_X);
            let torque = -SPIN_DRAG * (omega + spin);
            let mut calls = 0;
            let mut last = None;
            let result = joint::step_with_coupled_bushings_and_forcing(
                &model,
                &start,
                loads(point),
                &[closure(point)],
                &rows,
                &[coupled()],
                h,
                config(),
                |elapsed, mid, wrenches, effort| {
                    close(elapsed, half);
                    calls += 1;
                    let power = environment(&model, elapsed, mid, wrenches, effort)?;
                    last = Some((mid.clone(), power, wrenches[2]));
                    Ok(())
                },
                || false,
            )
            .unwrap();
            let step = &result.step.bilateral;
            let mid = &result.midpoint_state;
            let (last_state, power, wrench) = last.unwrap();
            assert_eq!(calls, step.residual_evaluations);
            assert_eq!(last_state, *mid);
            close(mid.base.twist_body.linear.x, u);
            close(mid.base.twist_body.angular.x, omega);
            for i in 0..2 {
                close(mid.joint_position[i], qm[i]);
                close(mid.joint_velocity[i], expected[i + 2]);
                close(
                    step.state_after.joint_position[i],
                    Q0[i] + h * expected[i + 2],
                );
                close(
                    step.state_after.joint_velocity[i],
                    2.0 * expected[i + 2] - start.joint_velocity[i],
                );
            }
            close(
                step.state_after.joint_velocity[2],
                step.state_after.joint_velocity[0],
            );
            close(
                step.state_after.base.twist_body.linear.x,
                2.0 * u - start.base.twist_body.linear.x,
            );
            close(
                step.state_after.base.twist_body.angular.x,
                2.0 * omega - start.base.twist_body.angular.x,
            );
            close(
                step.state_after.base.world_from_base.translation().x,
                2.0 + h * u,
            );
            close(wrench.force.x, force);
            close(wrench.torque.x, torque);
            close(
                step.midpoint.motion.external_linear_impulse_world.x,
                h * force,
            );
            let moment = step.midpoint.motion.external_angular_impulse_world;
            close(moment.x, h * torque);
            close(moment.y, 0.4 * h * force);
            close(moment.z, 0.7 * h * force);
            close(step.midpoint.motion.external_work, h * power);
            close(power, force * (u + vx) + torque * (omega + spin));
            let motor: [f64; 2] =
                core::array::from_fn(|i| HELD[i] - 2.0 * qm[i] - 0.3 * expected[i + 2]);
            close(
                step.midpoint.motion.actuator_work,
                h * (motor[0] * vx + motor[1] * spin),
            );
            let rates = [vx, spin];
            let loss = h
                * (0..2)
                    .map(|i| (0..2).map(|j| rates[i] * C[i][j] * rates[j]).sum::<f64>())
                    .sum::<f64>();
            close(step.midpoint.dissipated_energy, loss);
            let follower_mass = if point { 3.0 } else { 1.0 };
            let reaction = half
                * follower_mass
                * ((u - start.base.twist_body.linear.x) + (vx - start.joint_velocity[0]))
                / half;
            assert!(reaction.abs() > 1e-3);
            close(step.initial_impulse_ns[0], reaction);
            close(step.final_impulse_ns[0], reaction);
            if locked {
                let bushing_torque = K[1][0] * qm[0] + K[1][1] * qm[1] + C[1][0] * vx;
                let joint_reaction = half
                    * (2.0 * (omega - start.base.twist_body.angular.x) / half - motor[1]
                        + bushing_torque
                        - torque);
                close(result.step.initial_joint_multipliers[0], joint_reaction);
                close(result.step.final_joint_multipliers[0], joint_reaction);
                close(result.step.joint_position_residuals_after[0], 0.0);
                close(result.step.joint_velocity_residuals_after[0], 0.0);
            }
            assert!(step.initial_reaction.impulse_work > 1e-6);
            assert!(step.final_reaction.impulse_work < -1e-6);
            close(step.constraint_impulse_work, 0.0);
            close(step.mechanical_work_residual, 0.0);
            close(step.linear_momentum_balance_residual_world.norm(), 0.0);
            close(step.angular_momentum_balance_residual_world.norm(), 0.0);
            close(result.com_velocity_residual_world.norm(), 0.0);
            for residual in step
                .position_residual_after_m
                .iter()
                .chain(&step.velocity_residual_after_mps)
            {
                close(*residual, 0.0);
            }
            for point in &step.point_residual_after {
                close(point.relative_velocity_world_mps.norm(), 0.0);
            }
            assert_eq!(start.base, saved.base);
            assert_eq!(start.joint_position, saved.joint_position);
        }
    }
}

#[test]
fn forced_loop_diagonal_full_and_empty_joint_paths_are_exact() {
    let point = true;
    let model = model(point);
    let start = initial(point);
    let diagonal = [diagonal()];
    let full = [CoupledSpatialBushing::try_from_diagonal(diagonal[0], scales()).unwrap()];
    let forcing = |elapsed, stage: &FreeFloatingMotionState, w: &mut [Wrench], u: &mut [f64]| {
        environment(&model, elapsed, stage, w, u).map(|_| ())
    };
    let a = step_free_floating_bilateral_loops_with_forcing(
        &model,
        &start,
        loads(point),
        &[closure(point)],
        &diagonal,
        0.1,
        config(),
        forcing,
        || false,
    )
    .unwrap();
    let b = step_free_floating_bilateral_loops_with_coupled_bushings_and_forcing(
        &model,
        &start,
        loads(point),
        &[closure(point)],
        &full,
        0.1,
        config(),
        forcing,
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    let mixed = joint::step_with_bushings_and_forcing(
        &model,
        &start,
        loads(point),
        &[closure(point)],
        &[],
        &diagonal,
        0.1,
        config(),
        forcing,
        || false,
    )
    .unwrap();
    assert_eq!(mixed.step.bilateral, a.step);
    assert_eq!(mixed.midpoint_state, a.midpoint_state);
    assert_eq!(
        mixed.com_velocity_residual_world,
        a.com_velocity_residual_world
    );
    let no_bushing = step_free_floating_bilateral_loops_with_forcing(
        &model,
        &start,
        loads(point),
        &[closure(point)],
        &[],
        0.1,
        config(),
        forcing,
        || false,
    )
    .unwrap();
    let no_full = step_free_floating_bilateral_loops_with_coupled_bushings_and_forcing(
        &model,
        &start,
        loads(point),
        &[closure(point)],
        &[],
        0.1,
        config(),
        forcing,
        || false,
    )
    .unwrap();
    assert_eq!(no_bushing, no_full);
}

#[test]
fn forced_loop_trial_and_actual_impulse_verification_refuse_atomically() {
    let model = model(true);
    let start = initial(true);
    let saved = start.clone();
    let run = || {
        step_free_floating_bilateral_loops_with_coupled_bushings_and_forcing(
            &model,
            &start,
            loads(true),
            &[closure(true)],
            &[coupled()],
            0.2,
            config(),
            |elapsed, stage, w, u| environment(&model, elapsed, stage, w, u).map(|_| ()),
            || false,
        )
        .unwrap()
    };
    let reference = run();
    let last = reference.step.residual_evaluations;
    assert!(last > 2);
    for mode in 0..6 {
        let calls = Cell::new(0);
        let mut controls = config();
        if mode == 5 {
            controls.nonlinear.maximum_evaluations = last - 1;
        }
        let error = step_free_floating_bilateral_loops_with_coupled_bushings_and_forcing(
            &model,
            &start,
            loads(true),
            &[closure(true)],
            &[coupled()],
            0.2,
            controls,
            |elapsed, stage, w, u| {
                calls.set(calls.get() + 1);
                environment(&model, elapsed, stage, w, u)?;
                if (mode == 0 && calls.get() == 2) || (mode == 1 && calls.get() == last) {
                    return Err(FreeFloatingMotionError::InvalidLoad("environment domain"));
                }
                if mode == 2 && calls.get() == last {
                    w[2].force.x = f64::NAN;
                }
                if mode == 3 && calls.get() == last {
                    u[0] = f64::INFINITY;
                }
                Ok(())
            },
            || mode == 4 && calls.get() == last,
        )
        .unwrap_err();
        if mode == 5 {
            assert!(matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::EvaluationBudget { .. })
            ));
        } else if mode == 4 {
            assert!(matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Cancelled
                ))
            ));
        } else {
            assert!(matches!(
                error,
                BilateralLoopError::Motion(ImplicitSpatialError::Motion(_))
            ));
        }
        assert_eq!(start, saved);
    }
    assert_eq!(run(), reference);
}
