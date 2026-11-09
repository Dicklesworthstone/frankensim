//! G1: held actuation and damping share one continuous physical trajectory.

use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::driven::{
    DrivenViscousConnectorError, PreparedDrivenViscousConnector,
};

fn axial() -> [RigidBodyState; 2] {
    [
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::new(4.0, 0.0, 0.0), Vec3::ZERO),
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::new(-3.0, 0.0, 0.0), Vec3::ZERO),
    ]
}

#[test]
fn held_actuation_reaches_terminal_slip_and_integrates_work_and_heat() {
    let properties = scalar_properties();
    let states = axial();
    let mobility = 1.0 / 2.0 + 1.0 / 3.0;
    let coefficient = 3.0;
    let force = 6.0;
    let rate = coefficient * mobility;
    let terminal = force / coefficient;
    let difference = -3.0 - terminal;
    for duration in [0.0, 1e-6, 0.2, 10.0, 1e4] {
        let mut prepared = PreparedDrivenViscousConnector::new(
            states.map(|s| s.pose()), properties, [Vec3::ZERO; 2],
            diagonal([coefficient, 0.0, 0.0, 0.0, 0.0, 0.0]), duration, || false,
        ).unwrap();
        let result = prepared.advance(states, [force, 0.0, 0.0, 0.0, 0.0, 0.0], controls(2e-9), || false).unwrap();
        let slip = terminal + difference * f64::exp(-rate * duration);
        let integral = -f64::exp_m1(-rate * duration) / rate;
        let square_integral = -f64::exp_m1(-2.0 * rate * duration) / (2.0 * rate);
        let distance = terminal * duration + difference * integral;
        let heat = coefficient * (terminal * terminal * duration
            + 2.0 * terminal * difference * integral + difference * difference * square_integral);
        let impulse = (slip + 3.0) / mobility;
        close(result.relative_twist_after[0], slip, 2e-11);
        close(result.actuator_work_j, force * distance, 2e-9);
        close(result.dissipated_energy_j, heat, 2e-9);
        close(result.kinetic_balance_residual_j, 0.0, 2e-9);
        close(result.connector_balance_residual_j, 0.0, 2e-9);
        close(result.impulse.bodies[0].state_after.linear_momentum_world().x, 4.0 - impulse, 2e-11);
        close(result.impulse.bodies[1].state_after.linear_momentum_world().x, -3.0 + impulse, 2e-11);
        for i in 0..2 {
            assert_eq!(result.impulse.bodies[i].state_after.pose(), states[i].pose());
        }
        if duration > 0.0 && duration <= 0.2 {
            assert!(result.actuator_work_j < 0.0, "the actuator initially brakes the opposing slip");
            assert!(result.dissipated_energy_j > 0.0);
        }
        if duration > 10.0 {
            close(result.relative_twist_after[0], terminal, 2e-11);
        }
    }
}

#[test]
fn driven_undamped_directions_keep_their_physical_acceleration() {
    let properties = scalar_properties();
    let states = [body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO); 2];
    let duration = 2.0;
    let effort = [3.0, -2.0, 0.0, 0.0, 0.0, 0.8];
    for coefficient in [0.0, 3.0] {
        let mut prepared = PreparedDrivenViscousConnector::new(
            states.map(|s| s.pose()), properties, [Vec3::ZERO; 2],
            diagonal([coefficient, 0.0, 0.0, 0.0, 0.0, 0.0]), duration, || false,
        ).unwrap();
        let result = prepared.advance(states, effort, controls(2e-10), || false).unwrap();
        close(result.relative_twist_after[1], -2.0 * (1.0 / 2.0 + 1.0 / 3.0) * duration, 2e-13);
        close(result.relative_twist_after[5], 0.8 * (1.0 / 4.0 + 1.0 / 4.0) * duration, 2e-13);
        vector(result.impulse.linear_momentum_change_world_ns, Vec3::ZERO, 2e-13);
        vector(result.impulse.angular_momentum_change_world_nms, Vec3::ZERO, 2e-13);
        if coefficient == 0.0 {
            assert_eq!(result.dissipated_energy_j, 0.0);
            vector(result.impulse.bodies[1].linear_impulse_world, Vec3::new(6.0, -4.0, 0.0), 2e-13);
            vector(result.impulse.bodies[1].angular_impulse_body, Vec3::new(0.0, 0.0, 1.6), 2e-13);
        }
    }
}

fn forced_derivative(
    y: &[f64; 14], poses: [Pose; 2], arms: [Vec3; 2], c: [[f64; 6]; 6], effort: [f64; 6],
) -> [f64; 14] {
    let native: [f64; 13] = core::array::from_fn(|i| y[i]);
    let passive = derivative(&native, poses, arms, c);
    let force = poses[0].orientation().rotate_body_to_world(Vec3::new(effort[0], effort[1], effort[2]));
    let couple = poses[0].orientation().rotate_body_to_world(Vec3::new(effort[3], effort[4], effort[5]));
    let arm_b = poses[1].orientation().rotate_body_to_world(arms[1]);
    let displacement = poses[1].position_world().sub(poses[0].position_world());
    let torques = [
        poses[0].orientation().rotate_world_to_body(couple.add(displacement.add(arm_b).cross(force)).scale(-1.0)),
        poses[1].orientation().rotate_world_to_body(couple.add(arm_b.cross(force))),
    ];
    let mut drive = [0.0; 13];
    put(&mut drive, 0, force.scale(-1.0));
    put(&mut drive, 3, torques[0]);
    put(&mut drive, 6, force);
    put(&mut drive, 9, torques[1]);
    let mut result: [_; 14] = core::array::from_fn(|i| if i < 13 { passive[i] + drive[i] } else { 0.0 });
    // Independent actuator power in the original twelve momentum coordinates.
    result[13] = force.scale(-1.0).dot(get(&native, 0).scale(0.5))
        + force.dot(get(&native, 6).scale(1.0 / 3.0))
        + torques[0].dot(omega(get(&native, 3)))
        + torques[1].dot(omega(get(&native, 9)));
    result
}

#[test]
fn full_spatial_actuation_matches_continuous_newton_euler_and_subdivision() {
    let (states, properties, arms) = full_case();
    let poses = states.map(|s| s.pose());
    let c = full_matrix();
    let effort = [0.4, -0.7, 0.2, 0.8, 0.3, -0.5];
    let duration = 0.7;
    let mut reference = [0.0; 14];
    for i in 0..2 {
        let p = states[i].linear_momentum_world();
        let h = states[i].angular_momentum_body();
        reference[6*i..6*i+6].copy_from_slice(&[p.x, p.y, p.z, h.x, h.y, h.z]);
    }
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = forced_derivative(&reference, poses, arms, c, effort);
        let k2 = forced_derivative(&core::array::from_fn(|i| reference[i] + 0.5*h*k1[i]), poses, arms, c, effort);
        let k3 = forced_derivative(&core::array::from_fn(|i| reference[i] + 0.5*h*k2[i]), poses, arms, c, effort);
        let k4 = forced_derivative(&core::array::from_fn(|i| reference[i] + h*k3[i]), poses, arms, c, effort);
        for i in 0..14 {
            reference[i] += h * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        }
    }
    for steps in [1, 4] {
        let mut prepared = PreparedDrivenViscousConnector::new(poses, properties, arms, damping(c), duration / steps as f64, || false).unwrap();
        let mut current = states;
        let mut heat = 0.0;
        let mut work = 0.0;
        for _ in 0..steps {
            let result = prepared.advance(current, effort, controls(2e-9), || false).unwrap();
            current = result.impulse.bodies.map(|body| body.state_after);
            heat += result.dissipated_energy_j;
            work += result.actuator_work_j;
        }
        for i in 0..2 {
            vector(current[i].linear_momentum_world(), Vec3::new(reference[6*i], reference[6*i+1], reference[6*i+2]), 3e-11);
            vector(current[i].angular_momentum_body(), Vec3::new(reference[6*i+3], reference[6*i+4], reference[6*i+5]), 3e-11);
        }
        close(heat, reference[12], 3e-11);
        close(work, reference[13], 3e-11);
        vector(angular(current, Vec3::new(2.0, -3.0, 1.0)), angular(states, Vec3::new(2.0, -3.0, 1.0)), 2e-11);
    }
}

#[test]
fn unforced_limit_refusals_cancellation_and_reuse_keep_native_state() {
    let (states, properties, arms) = full_case();
    let poses = states.map(|s| s.pose());
    let law = damping(full_matrix());
    let passive = PreparedViscousConnector::new(poses, properties, arms, law.clone(), 0.3, || false).unwrap()
        .relax(states, controls(1e-9), || false).unwrap();
    let mut driven = PreparedDrivenViscousConnector::new(poses, properties, arms, law, 0.3, || false).unwrap();
    let result = driven.advance(states, [0.0; 6], controls(1e-9), || false).unwrap();
    for i in 0..2 {
        vector(result.impulse.bodies[i].state_after.linear_momentum_world(), passive.impulse.bodies[i].state_after.linear_momentum_world(), 2e-11);
        vector(result.impulse.bodies[i].state_after.angular_momentum_body(), passive.impulse.bodies[i].state_after.angular_momentum_body(), 2e-11);
    }
    close(result.dissipated_energy_j, passive.dissipated_energy_j, 2e-11);
    assert_eq!(result.actuator_work_j, 0.0);
    assert!(driven.advance(states, [f64::NAN; 6], controls(1e-9), || false).is_err());
    let mut changed = states;
    changed[0] = body(Vec3::ZERO, states[0].pose().orientation(), states[0].linear_momentum_world(), states[0].angular_momentum_body());
    assert!(matches!(driven.advance(changed, [0.0; 6], controls(1e-9), || false), Err(DrivenViscousConnectorError::Source(ViscousConnectorError::Invalid(_)))));
    for cancel_at in [0, 1, 5] {
        let mut polls = 0;
        assert!(driven.advance(states, [0.0; 6], controls(1e-9), || {
            let stop = polls == cancel_at;
            polls += 1;
            stop
        }).is_err());
    }
    assert_eq!(driven.advance(states, [0.0; 6], controls(1e-9), || false).unwrap(), result);
}
