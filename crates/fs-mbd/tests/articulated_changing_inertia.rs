//! Analytic changing-inertia trajectories for the original free-floating owner.
//! Both bodies share a COM and rotate about z, so individual spin momenta give
//! an independent scalar solution without a production mass matrix or solver.

use fs_ga::{Mat3, Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::step::{
    InertiaMotionError, advance, advance_with_effort,
};
use fs_mbd::articulated::motion::reconfigure::{ReconfigurationControls, ReconfigurationError};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use std::cell::Cell;

fn z(value: f64) -> Vec3 {
    Vec3::new(0., 0., value)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * (1. + expected.abs()),
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
fn model(carrier: f64, rotor: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., carrier),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(z(1.), None).unwrap(),
            sphere(1., rotor),
        ),
    ])
    .unwrap()
}
fn initial(omega: f64, joint_rate: f64) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(Se3::identity(), Twist::new(z(omega), Vec3::default())),
        joint_position: vec![0.2],
        joint_velocity: vec![joint_rate],
    }
}
fn controls(work: f64) -> ReconfigurationControls<'static> {
    ReconfigurationControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-11,
        linear_momentum_tolerance_ns: 2e-11,
        joint_momentum_tolerances: &[2e-11],
        maximum_parameter_work_j: work,
    }
}
fn kinetic(state: &FreeFloatingMotionState, carrier: f64, rotor: f64) -> f64 {
    let omega = state.base.twist_body.angular.z;
    0.5 * carrier * omega.powi(2) + 0.5 * rotor * (omega + state.joint_velocity[0]).powi(2)
}
fn momenta(state: &FreeFloatingMotionState, carrier: f64, rotor: f64) -> (f64, f64) {
    let omega = state.base.twist_body.angular.z;
    let joint = rotor * (omega + state.joint_velocity[0]);
    (carrier * omega + joint, joint)
}
fn angle(state: &FreeFloatingMotionState) -> f64 {
    let value = state.base.world_from_base.rotation().log().angular;
    close(value.x, 0.);
    close(value.y, 0.);
    close(state.base.twist_body.angular.x, 0.);
    close(state.base.twist_body.angular.y, 0.);
    close(state.base.twist_body.linear.norm(), 0.);
    close(state.base.world_from_base.translation().norm(), 0.);
    value.z
}
fn second_order(errors: [f64; 3], name: &str) {
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.5..4.5).contains(&ratio),
            "{name}: halving-step errors {errors:?} give ratio {ratio}"
        );
    }
}

// Integrals of f(t)=(p+k*t)/(j+a*t) and f(t)^2, from direct partial fractions.
fn integral(p: f64, k: f64, j: f64, a: f64, t: f64) -> f64 {
    k / a * t + (p - k * j / a) / a * ((j + a * t) / j).ln()
}
fn integral_square(p: f64, k: f64, j: f64, a: f64, t: f64) -> f64 {
    let b = p - k * j / a;
    (k / a).powi(2) * t
        + 2. * k * b / a.powi(2) * ((j + a * t) / j).ln()
        + b.powi(2) / a * (1. / j - 1. / (j + a * t))
}

#[test]
fn driven_changing_inertia_preserves_impulses_and_refines_angles_and_parameter_work() {
    let external_torque = 0.15;
    let joint_effort = 0.3;
    let exact_body_angle = integral(0.8, external_torque - joint_effort, 2., 1.2, 1.);
    let exact_joint_angle = integral(1.6, joint_effort, 1., 0.7, 1.) - exact_body_angle;
    let exact_parameter_work = -0.5
        * (1.2 * integral_square(0.8, external_torque - joint_effort, 2., 1.2, 1.)
            + 0.7 * integral_square(1.6, joint_effort, 1., 0.7, 1.));
    let initial_state = initial(0.4, 1.2);
    let exact_omega = (2.4 + external_torque - 1.6 - joint_effort) / 3.2;
    let exact_rate = (1.6 + joint_effort) / 1.7 - exact_omega;
    let exact_final_energy =
        0.5 * 3.2 * exact_omega.powi(2) + 0.5 * 1.7 * (exact_omega + exact_rate).powi(2);
    close(
        exact_final_energy
            - kinetic(&initial_state, 2., 1.)
            - exact_parameter_work
            - external_torque * exact_body_angle
            - joint_effort * exact_joint_angle,
        0.,
    );

    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = 1. / f64::from(n);
        let mut state = initial_state.clone();
        let (mut parameter_work, mut actuator_work, mut external_work, mut residual) =
            (0., 0., 0., 0.);
        for index in 0..n {
            let t = f64::from(index) * h;
            let previous = model(2. + 1.2 * t, 1. + 0.7 * t);
            let midpoint = model(2. + 1.2 * (t + h / 2.), 1. + 0.7 * (t + h / 2.));
            let next = model(2. + 1.2 * (t + h), 1. + 0.7 * (t + h));
            let wrenches = [
                Wrench::new(z(external_torque), Vec3::default()),
                Wrench::default(),
            ];
            let step = advance(
                &previous,
                &midpoint,
                &next,
                &state,
                FreeFloatingMotionLoads {
                    generalized_force: &[joint_effort],
                    external_body_wrench: &wrenches,
                    gravity_world: Vec3::default(),
                },
                h,
                controls(10.),
                || false,
            )
            .unwrap();
            let end = step.state_after();
            let (pb, pq) = momenta(end, 2. + 1.2 * (t + h), 1. + 0.7 * (t + h));
            close(pb, 2.4 + external_torque * (t + h));
            close(pq, 1.6 + joint_effort * (t + h));
            close(
                step.motion.external_angular_impulse_world.z,
                external_torque * h,
            );
            close(step.linear_momentum_residual.norm(), 0.);
            close(step.angular_momentum_residual.norm(), 0.);
            close(
                step.motion.actuator_work,
                joint_effort * (end.joint_position[0] - state.joint_position[0]),
            );
            close(
                step.motion.external_work,
                external_torque * (angle(end) - angle(&state)),
            );
            let change = kinetic(end, 2. + 1.2 * (t + h), 1. + 0.7 * (t + h))
                - kinetic(&state, 2. + 1.2 * t, 1. + 0.7 * t);
            close(
                step.kinetic_work_residual,
                change
                    - step.net_parameter_work_j
                    - step.motion.actuator_work
                    - step.motion.external_work,
            );
            parameter_work += step.net_parameter_work_j;
            actuator_work += step.motion.actuator_work;
            external_work += step.motion.external_work;
            residual += step.kinetic_work_residual;
            state = step.into_state_after();
        }
        close(state.base.twist_body.angular.z, exact_omega);
        close(state.joint_velocity[0], exact_rate);
        close(
            kinetic(&state, 3.2, 1.7)
                - kinetic(&initial_state, 2., 1.)
                - parameter_work
                - actuator_work
                - external_work,
            residual,
        );
        errors.push([
            (angle(&state) - exact_body_angle).abs(),
            (state.joint_position[0] - initial_state.joint_position[0] - exact_joint_angle).abs(),
            (parameter_work - exact_parameter_work).abs(),
        ]);
    }
    for (i, name) in [
        "driven base angle",
        "driven joint angle",
        "driven parameter work",
    ]
    .iter()
    .enumerate()
    {
        second_order([errors[0][i], errors[1][i], errors[2][i]], name);
        assert!(errors[2][i] < 2e-5);
    }
}

#[test]
fn identity_gross_work_and_pre_callback_refusals_keep_native_semantics() {
    let old = model(2., 1.);
    let state = initial(0.4, 1.2);
    let retained = (old.clone(), state.clone());
    let wrenches = [Wrench::new(z(0.15), Vec3::default()), Wrench::default()];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.3],
        external_body_wrench: &wrenches,
        gravity_world: Vec3::default(),
    };
    let original = step_free_floating(&old, &state, loads, 0.01, || false).unwrap();
    let identity = advance(&old, &old, &old, &state, loads, 0.01, controls(0.), || {
        false
    })
    .unwrap();
    assert_eq!(identity.motion, original);
    assert_eq!(identity.state_after(), &original.state_after);
    assert_eq!(identity.net_parameter_work_j, 0.);
    assert_eq!(identity.absolute_parameter_work_j, 0.);

    let midpoint = model(4., 3.);
    let zero_wrenches = [Wrench::default(); 2];
    let unloaded = FreeFloatingMotionLoads {
        generalized_force: &[0.],
        external_body_wrench: &zero_wrenches,
        gravity_world: Vec3::default(),
    };
    assert!(matches!(
        advance(
            &old,
            &midpoint,
            &old,
            &state,
            unloaded,
            0.01,
            controls(1.),
            || false
        ),
        Err(InertiaMotionError::Reconfiguration(
            ReconfigurationError::ParameterWork { .. }
        ))
    ));
    let cycle = advance(
        &old,
        &midpoint,
        &old,
        &state,
        unloaded,
        0.01,
        controls(2.),
        || false,
    )
    .unwrap();
    close(cycle.first_transfer.parameter_work_j, -14. / 15.);
    close(cycle.final_transfer.parameter_work_j, 14. / 15.);
    close(cycle.net_parameter_work_j, 0.);
    close(cycle.absolute_parameter_work_j, 28. / 15.);
    close(cycle.state_after().base.twist_body.angular.z, 0.4);
    close(cycle.state_after().joint_velocity[0], 1.2);
    close(angle(cycle.state_after()), 0.002);
    close(
        cycle.state_after().joint_position[0],
        state.joint_position[0] + 0.01 / 3.,
    );

    let changed_mass = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2.1, 2.),
        ),
        old.links()[1].clone(),
    ])
    .unwrap();
    for (next, duration) in [(&changed_mass, 0.01), (&old, 0.), (&old, f64::NAN)] {
        let mut callbacks = 0;
        let refused = advance_with_effort(
            &old,
            &old,
            next,
            &state,
            loads,
            duration,
            controls(10.),
            |_, _, _| {
                callbacks += 1;
                Ok(())
            },
            || false,
        );
        assert!(refused.is_err());
        assert_eq!(
            callbacks, 0,
            "inventory and duration admission must precede effort callbacks"
        );
    }
    let cancelled = Cell::new(false);
    let mut stages = Vec::new();
    let refused = advance_with_effort(
        &old,
        &midpoint,
        &old,
        &state,
        loads,
        0.01,
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
        Err(InertiaMotionError::Motion(
            FreeFloatingMotionError::Cancelled
        ))
    ));
    assert_eq!(stages, vec![0., 0.005]);
    assert_eq!((old, state), retained);
}

#[test]
fn damping_uses_virtual_midpoint_inertia_rates_and_refines_physical_work() {
    // pb=0, pq=1/(1+t), Jb=Jr=1+t. Hence omega=-1/(1+t)^2,
    // qdot=2/(1+t)^2, K=1/(1+t)^3. The two physical work channels
    // integrate to Wparameter=-7/24 and Wactuator=-7/12 at t=1.
    let original = initial(-1., 2.);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = 1. / f64::from(n);
        let mut state = original.clone();
        let (mut parameter_work, mut actuator_work, mut residual) = (0., 0., 0.);
        for index in 0..n {
            let t = f64::from(index) * h;
            let j0 = 1. + t;
            let jm = j0 + h / 2.;
            let j1 = j0 + h;
            let previous = model(j0, j0);
            let midpoint = model(jm, jm);
            let next = model(j1, j1);
            let (_, old_pq) = momenta(&state, j0, j0);
            let mut stages = Vec::new();
            let step = advance_with_effort(
                &previous,
                &midpoint,
                &next,
                &state,
                FreeFloatingMotionLoads {
                    generalized_force: &[0.],
                    external_body_wrench: &[Wrench::default(); 2],
                    gravity_world: Vec3::default(),
                },
                h,
                controls(10.),
                |elapsed, stage, efforts| {
                    assert_eq!(elapsed, if stages.is_empty() { 0. } else { h / 2. });
                    assert_eq!(efforts, &[0.]);
                    if stages.is_empty() {
                        close(stage.joint_velocity[0], 2. * old_pq / jm);
                    }
                    stages.push((elapsed, stage.joint_velocity[0]));
                    efforts[0] = -0.5 * stage.joint_velocity[0];
                    Ok(())
                },
                || false,
            )
            .unwrap();
            assert_eq!(stages.len(), 2);
            close(step.motion.actuator_work, -0.5 * h * stages[1].1.powi(2));
            close(step.motion.external_work, 0.);
            close(step.motion.internal_body_work, 0.);
            close(step.motion.external_angular_impulse_world.norm(), 0.);
            close(step.linear_momentum_residual.norm(), 0.);
            close(step.angular_momentum_residual.norm(), 0.);
            let end = step.state_after();
            close(momenta(end, j1, j1).0, 0.);
            let independent_residual = kinetic(end, j1, j1)
                - kinetic(&state, j0, j0)
                - step.net_parameter_work_j
                - step.motion.actuator_work;
            close(step.kinetic_work_residual, independent_residual);
            close(
                step.kinetic_work_residual,
                step.motion.kinetic_work_residual,
            );
            parameter_work += step.net_parameter_work_j;
            actuator_work += step.motion.actuator_work;
            residual += independent_residual;
            state = step.into_state_after();
        }
        let end_energy = kinetic(&state, 2., 2.);
        close(end_energy - 1. - parameter_work - actuator_work, residual);
        assert!(
            residual < 0.,
            "native explicit midpoint has a resolved numerical work defect here"
        );
        errors.push([
            (state.base.twist_body.angular.z + 0.25).abs(),
            (state.joint_velocity[0] - 0.5).abs(),
            (angle(&state) + 0.5).abs(),
            (state.joint_position[0] - original.joint_position[0] - 1.).abs(),
            (parameter_work + 7. / 24.).abs(),
            (actuator_work + 7. / 12.).abs(),
            (end_energy - 0.125).abs(),
            residual.abs(),
        ]);
    }
    let labels = [
        "damped base rate",
        "damped joint rate",
        "damped base angle",
        "damped joint angle",
        "damped parameter work",
        "damped actuator work",
        "damped kinetic energy",
        "numerical work residual",
    ];
    for (i, name) in labels.iter().enumerate() {
        second_order([errors[0][i], errors[1][i], errors[2][i]], name);
        assert!(errors[2][i] < 3e-4);
    }
}
