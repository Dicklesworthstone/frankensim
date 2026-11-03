//! Analytic gear-train motion with changing inertia and actual stage reactions.
//! A carrier and two coaxial rotors obey 2*q1+q2=0. Independent two-coordinate
//! momentum equations give trajectories, reaction forces and physical work.

use fs_ga::{Mat3, Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::constraints::step::{
    ConstrainedInertiaMotionError, advance, advance_with_effort,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintControls, ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::reconfigure::{ReconfigurationControls, step as free_step};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};
use std::cell::Cell;

fn z(value: f64) -> Vec3 {
    Vec3::new(0., 0., value)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * (1. + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn sphere(mass: f64, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [moment, 0., 0., 0., moment, 0., 0., 0., moment],
        },
    )
    .unwrap()
}
fn model(carrier: f64, first: f64, second: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., carrier),
        ),
        Link::new(
            "first rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(z(1.), None).unwrap(),
            sphere(1., first),
        ),
        Link::new(
            "second rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(z(1.), None).unwrap(),
            sphere(1., second),
        ),
    ])
    .unwrap()
}
fn state(omega: f64, rate: f64) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(Se3::identity(), Twist::new(z(omega), Vec3::default())),
        joint_position: vec![0.2, -0.4],
        joint_velocity: vec![rate, -2. * rate],
    }
}
fn gear() -> JointConstraint {
    JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![2., 1.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }
}
fn controls(work: f64) -> ConstraintControls<'static> {
    ConstraintControls {
        reconfiguration: ReconfigurationControls {
            relative_momentum_tolerance: 2e-12,
            angular_momentum_tolerance_nms: 1e-9,
            linear_momentum_tolerance_ns: 1e-9,
            joint_momentum_tolerances: &[1e-9, 1e-9],
            maximum_parameter_work_j: work,
        },
        maximum_reaction_work_j: 1e-8,
        minimum_reciprocal_condition: 1e-12,
    }
}
fn kinetic(s: &FreeFloatingMotionState, jb: f64, j1: f64, j2: f64) -> f64 {
    let w = s.base.twist_body.angular.z;
    0.5 * jb * w * w
        + 0.5 * j1 * (w + s.joint_velocity[0]).powi(2)
        + 0.5 * j2 * (w + s.joint_velocity[1]).powi(2)
}
fn momenta(s: &FreeFloatingMotionState, jb: f64, j1: f64, j2: f64) -> (f64, f64) {
    let w = s.base.twist_body.angular.z;
    let p1 = j1 * (w + s.joint_velocity[0]);
    let p2 = j2 * (w + s.joint_velocity[1]);
    (jb * w + p1 + p2, p1 - 2. * p2)
}
fn angle(s: &FreeFloatingMotionState) -> f64 {
    close(s.base.world_from_base.translation().norm(), 0.);
    close(s.base.twist_body.linear.norm(), 0.);
    s.base.world_from_base.rotation().log().angular.z
}
fn admit_gear(s: &FreeFloatingMotionState) {
    close(2. * s.joint_position[0] + s.joint_position[1], 0.);
    close(2. * s.joint_velocity[0] + s.joint_velocity[1], 0.);
}
fn integral(p: f64, k: f64, j: f64, a: f64) -> f64 {
    k / a + (p - k * j / a) / a * ((j + a) / j).ln()
}
fn integral_square(p: f64, k: f64, j: f64, a: f64) -> f64 {
    let b = p - k * j / a;
    (k / a).powi(2)
        + 2. * k * b / a.powi(2) * ((j + a) / j).ln()
        + b.powi(2) / a * (1. / j - 1. / (j + a))
}
fn second_order(errors: [f64; 3], label: &str) {
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.4..4.6).contains(&ratio),
            "{label}: {errors:?}, ratio {ratio}"
        );
    }
}

#[test]
fn driven_gear_reactions_and_nonuniform_inertia_refine_to_the_scalar_solution() {
    // Jb=2, J1=1+t, J2=3(1+t). The reduced mass is
    // [[2+4k,-5k],[-5k,13k]], k=1+t, in coordinates (base angle,q1).
    // pb=-3.6+0.15t, ps=13.6+0.7t; w=(13pb+5ps)/(26+27k),
    // s=ps/(13k)+(5/13)w. All expectations below follow these scalars.
    let original = state(0.4, 1.2);
    let exact_angle = integral(21.2, 5.45, 53., 27.);
    let exact_joint = integral(13.6, 0.7, 1., 1.) / 13. + 5. * exact_angle / 13.;
    let exact_parameter = -0.5 / 13.
        * (27. * integral_square(21.2, 5.45, 53., 27.) + integral_square(13.6, 0.7, 1., 1.));
    let rows = [gear()];
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = 1. / f64::from(n);
        let mut current = original.clone();
        let (mut parameter, mut actuator, mut external, mut reaction, mut residual) =
            (0., 0., 0., 0., 0.);
        let mut saw_transfer_reaction = false;
        for index in 0..n {
            let t = f64::from(index) * h;
            let k = 1. + t;
            let km = k + h / 2.;
            let previous = model(2., k, 3. * k);
            let midpoint = model(2., km, 3. * km);
            let next = model(2., k + h, 3. * (k + h));
            let loads = FreeFloatingMotionLoads {
                generalized_force: &[0.3, -0.2],
                external_body_wrench: &[
                    Wrench::new(z(0.15), Vec3::default()),
                    Wrench::default(),
                    Wrench::default(),
                ],
                gravity_world: Vec3::default(),
            };
            let step = advance(
                &previous,
                &midpoint,
                &next,
                &current,
                loads,
                h,
                &rows,
                controls(10.),
                || false,
            )
            .unwrap();
            let end = step.state_after();
            admit_gear(end);
            let (pb, ps) = momenta(end, 2., k + h, 3. * (k + h));
            close(pb, -3.6 + 0.15 * (t + h));
            close(ps, 13.6 + 0.7 * (t + h));
            close(
                end.base.twist_body.angular.z,
                (21.2 + 5.45 * (t + h)) / (53. + 27. * (t + h)),
            );
            close(
                end.joint_velocity[0],
                ps / (13. * (k + h)) + 5. * end.base.twist_body.angular.z / 13.,
            );
            assert_eq!(step.smooth_reactions.len(), 2);
            // Solve the independent scalar acceleration constraint. The frozen
            // midpoint model has no gyroscopic bias in this coaxial motion.
            let expected_lambda = (-3.2 + 0.45 * km) / (26. + 27. * km);
            for (stage, elapsed) in step.smooth_reactions.iter().zip([0., h / 2.]) {
                assert_eq!(stage.elapsed_s, elapsed);
                close(stage.multipliers[0], expected_lambda);
                close(stage.joint_effort[0], 2. * expected_lambda);
                close(stage.joint_effort[1], expected_lambda);
                close(stage.power_w, 0.);
            }
            saw_transfer_reaction |= step.first_transfer.multipliers[0].abs() > 1e-8;
            close(
                step.caller_actuator_work_j,
                0.7 * (end.joint_position[0] - current.joint_position[0]),
            );
            close(
                step.motion.external_work,
                0.15 * (angle(end) - angle(&current)),
            );
            close(step.net_reaction_work_j, 0.);
            close(step.linear_momentum_residual.norm(), 0.);
            close(step.angular_momentum_residual.norm(), 0.);
            let independent = kinetic(end, 2., k + h, 3. * (k + h))
                - kinetic(&current, 2., k, 3. * k)
                - step.net_parameter_work_j
                - step.caller_actuator_work_j
                - step.motion.external_work
                - step.net_reaction_work_j;
            close(step.kinetic_work_residual, independent);
            parameter += step.net_parameter_work_j;
            actuator += step.caller_actuator_work_j;
            external += step.motion.external_work;
            reaction += step.net_reaction_work_j;
            residual += step.kinetic_work_residual;
            current = step.into_state_after();
        }
        assert!(
            saw_transfer_reaction,
            "nonuniform inertia changes must redistribute joint momentum"
        );
        close(
            kinetic(&current, 2., 2., 6.)
                - kinetic(&original, 2., 1., 3.)
                - parameter
                - actuator
                - external
                - reaction,
            residual,
        );
        close(
            kinetic(&current, 2., 2., 6.)
                - kinetic(&original, 2., 1., 3.)
                - exact_parameter
                - 0.7 * exact_joint
                - 0.15 * exact_angle,
            0.,
        );
        errors.push([
            (angle(&current) - exact_angle).abs(),
            (current.joint_position[0] - original.joint_position[0] - exact_joint).abs(),
            (parameter - exact_parameter).abs(),
        ]);
    }
    for (i, label) in ["base angle", "joint angle", "parameter work"]
        .iter()
        .enumerate()
    {
        second_order([errors[0][i], errors[1][i], errors[2][i]], label);
        assert!(errors[2][i] < 1e-4);
    }
}

#[test]
fn damped_gear_uses_actual_stage_reactions_and_refines_work_and_motion() {
    // Scale all three inertias by k=1+t. With pb=0, ps0=53/6 and
    // tau1=-(53/6)*qdot1, ps=ps0/k, qdot1=1/k² and w=5/(6k²).
    // K=53/(12k³); at t=1, Wparameter=-371/288, Wactuator=-371/144.
    let original = state(5. / 6., 1.);
    let rows = [gear()];
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = 1. / f64::from(n);
        let mut current = original.clone();
        let (mut parameter, mut actuator, mut residual) = (0., 0., 0.);
        for index in 0..n {
            let k = 1. + f64::from(index) * h;
            let km = k + h / 2.;
            let previous = model(2. * k, k, 3. * k);
            let midpoint = model(2. * km, km, 3. * km);
            let next = model(2. * (k + h), k + h, 3. * (k + h));
            let (_, old_ps) = momenta(&current, 2. * k, k, 3. * k);
            let mut stages = Vec::new();
            let step = advance_with_effort(
                &previous,
                &midpoint,
                &next,
                &current,
                FreeFloatingMotionLoads {
                    generalized_force: &[0., 0.],
                    external_body_wrench: &[Wrench::default(); 3],
                    gravity_world: Vec3::default(),
                },
                h,
                &rows,
                controls(10.),
                |elapsed, stage, efforts| {
                    assert_eq!(elapsed, if stages.is_empty() { 0. } else { h / 2. });
                    assert_eq!(efforts, &[0., 0.]);
                    admit_gear(stage);
                    if stages.is_empty() {
                        close(stage.joint_velocity[0], 6. * old_ps / (53. * km));
                    }
                    stages.push(stage.joint_velocity[0]);
                    efforts[0] = -(53. / 6.) * stage.joint_velocity[0];
                    Ok(())
                },
                || false,
            )
            .unwrap();
            assert_eq!(stages.len(), 2);
            close(
                step.caller_actuator_work_j,
                -(53. / 6.) * h * stages[1].powi(2),
            );
            close(step.net_reaction_work_j, 0.);
            close(step.motion.external_work, 0.);
            close(step.motion.internal_body_work, 0.);
            close(step.linear_momentum_residual.norm(), 0.);
            close(step.angular_momentum_residual.norm(), 0.);
            let end = step.state_after();
            admit_gear(end);
            close(momenta(end, 2. * (k + h), k + h, 3. * (k + h)).0, 0.);
            let independent = kinetic(end, 2. * (k + h), k + h, 3. * (k + h))
                - kinetic(&current, 2. * k, k, 3. * k)
                - step.net_parameter_work_j
                - step.caller_actuator_work_j
                - step.net_reaction_work_j;
            close(step.kinetic_work_residual, independent);
            parameter += step.net_parameter_work_j;
            actuator += step.caller_actuator_work_j;
            residual += independent;
            current = step.into_state_after();
        }
        close(
            kinetic(&current, 4., 2., 6.) - 53. / 12. - parameter - actuator,
            residual,
        );
        assert!(
            residual < 0.,
            "native midpoint numerical defect stays visible"
        );
        errors.push([
            (current.joint_velocity[0] - 0.25).abs(),
            (current.base.twist_body.angular.z - 5. / 24.).abs(),
            (current.joint_position[0] - original.joint_position[0] - 0.5).abs(),
            (angle(&current) - 5. / 12.).abs(),
            (parameter + 371. / 288.).abs(),
            (actuator + 371. / 144.).abs(),
            residual.abs(),
        ]);
    }
    for (i, label) in [
        "joint rate",
        "base rate",
        "joint angle",
        "base angle",
        "parameter work",
        "actuator work",
        "numerical defect",
    ]
    .iter()
    .enumerate()
    {
        second_order([errors[0][i], errors[1][i], errors[2][i]], label);
        assert!(errors[2][i] < 2e-3);
    }
}

#[test]
fn empty_rows_gross_work_and_midpoint_cancellation_keep_atomic_semantics() {
    let old = model(2., 1., 3.);
    let midpoint = model(4., 2., 6.);
    let initial = state(0.4, 1.2);
    let retained = initial.clone();
    let rows = [gear()];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0., 0.],
        external_body_wrench: &[Wrench::default(); 3],
        gravity_world: Vec3::default(),
    };
    let free = free_step::advance(
        &old,
        &midpoint,
        &old,
        &initial,
        loads,
        0.01,
        controls(10.).reconfiguration,
        || false,
    )
    .unwrap();
    let empty = advance(
        &old,
        &midpoint,
        &old,
        &initial,
        loads,
        0.01,
        &[],
        controls(10.),
        || false,
    )
    .unwrap();
    assert_eq!(empty.motion, free.motion);
    assert_eq!(empty.state_after(), free.state_after());
    assert_eq!(empty.net_parameter_work_j, free.net_parameter_work_j);
    assert_eq!(
        empty.absolute_parameter_work_j,
        free.absolute_parameter_work_j
    );
    assert!(
        advance(
            &old,
            &midpoint,
            &old,
            &initial,
            loads,
            0.01,
            &rows,
            controls(5.),
            || false
        )
        .is_err()
    );
    let cycle = advance(
        &old,
        &midpoint,
        &old,
        &initial,
        loads,
        0.01,
        &rows,
        controls(8.),
        || false,
    )
    .unwrap();
    close(cycle.net_parameter_work_j, 0.);
    close(cycle.absolute_parameter_work_j, 7.44);
    close(cycle.state_after().base.twist_body.angular.z, 0.4);
    close(cycle.state_after().joint_velocity[0], 1.2);
    admit_gear(cycle.state_after());
    let cancelled = Cell::new(false);
    let mut stages = Vec::new();
    let refused = advance_with_effort(
        &old,
        &midpoint,
        &old,
        &initial,
        loads,
        0.01,
        &rows,
        controls(10.),
        |elapsed, _, _| {
            stages.push(elapsed);
            if elapsed > 0. {
                cancelled.set(true);
            }
            Ok(())
        },
        || cancelled.get(),
    );
    assert!(matches!(
        refused,
        Err(ConstrainedInertiaMotionError::Constraint(
            ConstraintError::Cancelled
        ))
    ));
    assert_eq!(stages, vec![0., 0.005]);
    assert_eq!(initial, retained);
}

fn locked_slider() -> (ArticulatedModel, FreeFloatingMotionState, JointConstraint) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(3., 2.),
        ),
        Link::new(
            "locked slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1., 0., 0.),
                Some(JointLimits::new(-1., 1., 0.1, 0.1).unwrap()),
            )
            .unwrap(),
            sphere(2., 1.),
        ),
    ])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(1., 0., 0.)),
        ),
        joint_position: vec![0.2],
        joint_velocity: vec![0.],
    };
    let row = JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![1.],
        position: 0.2,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    };
    (model, state, row)
}

fn lock_controls() -> ConstraintControls<'static> {
    ConstraintControls {
        reconfiguration: ReconfigurationControls {
            joint_momentum_tolerances: &[1e-9],
            ..controls(0.).reconfiguration
        },
        ..controls(0.)
    }
}

#[test]
fn passive_lock_transmits_external_force_without_consuming_motor_effort() {
    let (model, initial, row) = locked_slider();
    let saved = initial.clone();
    let force = Wrench::new(Vec3::default(), Vec3::new(100., 0., 0.));
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.],
        external_body_wrench: &[Wrench::default(), force],
        gravity_world: Vec3::default(),
    };
    let dt = 0.2;
    let mut stages = Vec::new();
    let result = advance_with_effort(
        &model,
        &model,
        &model,
        &initial,
        loads,
        dt,
        &[row.clone()],
        lock_controls(),
        |elapsed, stage, efforts| {
            stages.push(elapsed);
            assert_eq!(efforts, &[0.]);
            close(stage.joint_velocity[0], 0.);
            Ok(())
        },
        || false,
    )
    .unwrap();
    // The locked total mass is 5 kg. Its 100 N external load produces 20 m/s²
    // carrier acceleration and -60 N at the passive joint, despite its 0.1 N
    // actuator rating. The unconstrained 50 m/s² joint response is never an
    // admitted physical trajectory.
    assert_eq!(stages, vec![0., dt / 2.]);
    for reaction in &result.smooth_reactions {
        close(reaction.joint_effort[0], -60.);
        close(reaction.power_w, 0.);
        close(reaction.caller_actuator_power_w, 0.);
    }
    let endpoint = result.state_after();
    close(endpoint.base.twist_body.linear.x, 5.);
    close(endpoint.base.world_from_base.translation().x, 0.6);
    close(endpoint.joint_position[0], 0.2);
    close(endpoint.joint_velocity[0], 0.);
    close(result.motion.linear_momentum_change_world.x, 20.);
    close(result.motion.external_linear_impulse_world.x, 20.);
    close(result.motion.angular_momentum_change_world.norm(), 0.);
    close(result.motion.actuator_work, 0.);
    close(result.motion.constraint_work, 0.);
    close(result.motion.external_work, 60.);
    close(result.motion.kinetic_energy_change, 60.);
    close(result.motion.kinetic_work_residual, 0.);
    close(result.kinetic_work_residual, 0.);
    assert_eq!(initial, saved);
    let refusal = advance_with_effort(
        &model,
        &model,
        &model,
        &initial,
        loads,
        dt,
        &[row],
        lock_controls(),
        |_, _, efforts| {
            efforts[0] = 0.2;
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        refusal,
        ConstrainedInertiaMotionError::Motion(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointEffortLimitViolation { .. }
        ))
    ));
}

#[test]
fn finite_tangency_defect_retains_signed_work_separate_from_actuation() {
    let (model, mut initial, mut row) = locked_slider();
    initial.joint_velocity[0] = 0.02;
    row.position_tolerance = 0.01;
    row.velocity_tolerance = 0.1;
    let saved = initial.clone();
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.],
        external_body_wrench: &[
            Wrench::default(),
            Wrench::new(Vec3::default(), Vec3::new(100., 0., 0.)),
        ],
        gravity_world: Vec3::default(),
    };
    let controls = ConstraintControls {
        maximum_reaction_work_j: 0.3,
        ..lock_controls()
    };
    let result = advance(
        &model,
        &model,
        &model,
        &initial,
        loads,
        0.2,
        &[row.clone()],
        controls,
        || false,
    )
    .unwrap();
    // The admitted 0.02 m/s row defect has actual signed power -60*.02 W.
    // External work is 100*(3+.02)*.2 J; physical kinetic gain is 60.16 J.
    close(result.motion.actuator_work, 0.);
    close(result.caller_actuator_work_j, 0.);
    close(result.motion.constraint_work, -0.24);
    close(result.smooth_reaction_work_j, -0.24);
    close(result.absolute_reaction_work_j, 0.24);
    close(result.motion.external_work, 60.4);
    close(result.motion.kinetic_energy_change, 60.16);
    close(result.motion.kinetic_work_residual, 0.);
    close(result.kinetic_work_residual, 0.);
    close(result.state_after().joint_velocity[0], 0.02);
    close(result.state_after().joint_position[0], 0.204);
    let refusal = advance(
        &model,
        &model,
        &model,
        &initial,
        loads,
        0.2,
        &[row],
        ConstraintControls {
            maximum_reaction_work_j: 0.1,
            ..controls
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        refusal,
        ConstrainedInertiaMotionError::ReactionWorkBudget { .. }
    ));
    assert_eq!(initial, saved);
}
