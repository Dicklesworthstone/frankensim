//! Full passive translation/rotation coupling inside native closed mechanisms.
//!
//! G1 reduced two-mode midpoint oracles retain both carrier recoils, rigid
//! reactions, finite storage and loss. G0/G4 retain parity and atomic refusal.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, BilateralLoopStep, DistanceLoop,
    PointCoincidenceLoop, step_free_floating_bilateral_loops,
    step_free_floating_bilateral_loops_with_bushings,
    step_free_floating_bilateral_loops_with_coupled_bushings,
    step_free_floating_bilateral_loops_with_coupled_bushings_and_effort,
    step_free_floating_bilateral_loops_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
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
