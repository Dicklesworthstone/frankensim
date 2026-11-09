//! Physical full spatial damping, with independent continuous momentum oracles.

use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix, PreparedViscousConnector,
    ViscousConnectorControls, ViscousConnectorError,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

#[path = "internal_viscous_connector/driven.rs"]
mod driven;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn controls(tolerance: f64) -> ViscousConnectorControls {
    ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::ZERO,
            linear_tolerance_ns: tolerance,
            angular_tolerance_nms: tolerance,
            work_tolerance_j: tolerance,
        },
        linear_speed_tolerance_m_s: tolerance,
        angular_speed_tolerance_rad_s: tolerance,
    }
}
fn damping(values: [[f64; 6]; 6]) -> PassiveBushingMatrix {
    PassiveBushingMatrix::try_new(
        values,
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
    )
    .unwrap()
}
fn diagonal(values: [f64; 6]) -> PassiveBushingMatrix {
    damping(core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { values[i] } else { 0.0 })
    }))
}
fn body(position: Vec3, orientation: UnitQuaternion, p: Vec3, h: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::new(position, orientation).unwrap(), p, h).unwrap()
}
fn scalar_properties() -> [MassProperties; 2] {
    [2.0, 3.0].map(|mass| MassProperties::new(mass, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap())
}
fn angular(states: [RigidBodyState; 2], origin: Vec3) -> Vec3 {
    states.into_iter().fold(Vec3::ZERO, |sum, state| {
        sum.add(
            state
                .pose()
                .position_world()
                .sub(origin)
                .cross(state.linear_momentum_world()),
        )
        .add(
            state
                .pose()
                .orientation()
                .rotate_body_to_world(state.angular_momentum_body()),
        )
    })
}

#[test]
fn axial_limit_preserves_recoil_and_integrates_weak_and_stiff_decay() {
    let properties = scalar_properties();
    let states = [
        body(
            Vec3::new(-1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::ZERO,
        ),
        body(
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(-3.0, 0.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let mobility = 1.0 / 2.0 + 1.0 / 3.0;
    let coefficient = 3.0;
    let initial_slip = -3.0;
    for duration in [0.0, 1e-12, 0.2, 10.0, 1e4] {
        let prepared = PreparedViscousConnector::new(
            states.map(|s| s.pose()),
            properties,
            [Vec3::ZERO; 2],
            diagonal([coefficient, 0.0, 0.0, 0.0, 0.0, 0.0]),
            duration,
            || false,
        )
        .unwrap();
        let result = prepared.relax(states, controls(2e-11), || false).unwrap();
        let final_slip = initial_slip * f64::exp(-coefficient * mobility * duration);
        let impulse = (final_slip - initial_slip) / mobility;
        let heat = initial_slip * initial_slip / (2.0 * mobility)
            * -f64::exp_m1(-2.0 * coefficient * mobility * duration);
        close(
            result.impulse.bodies[0]
                .state_after
                .linear_momentum_world()
                .x,
            4.0 - impulse,
            2e-12,
        );
        close(
            result.impulse.bodies[1]
                .state_after
                .linear_momentum_world()
                .x,
            -3.0 + impulse,
            2e-12,
        );
        close(result.relative_twist_after[0], final_slip, 2e-12);
        assert!(result.relative_twist_after[0] <= 2e-14);
        close(result.dissipated_energy_j, heat, 2e-12);
        close(result.kinetic_heat_residual_j, 0.0, 2e-12);
        for i in 0..2 {
            assert_eq!(
                result.impulse.bodies[i].state_after.pose(),
                states[i].pose()
            );
            assert_eq!(
                result.impulse.bodies[i].state_after.angular_momentum_body(),
                Vec3::ZERO
            );
        }
    }
}

#[test]
fn transverse_shear_includes_transport_torque_and_rotational_mobility() {
    let properties = scalar_properties();
    let states = [
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO),
        body(
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(0.0, 9.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let duration = 0.7;
    let coefficient = 5.0;
    // B's force acts at x=2; transport to A supplies the opposite 2*J couple.
    // This spin recoil contributes 2^2/I_A,z to the same shear mobility.
    let mobility = 1.0 / 2.0 + 1.0 / 3.0 + 4.0 / 4.0;
    let final_slip = 3.0 * f64::exp(-coefficient * mobility * duration);
    let impulse = (final_slip - 3.0) / mobility;
    let prepared = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        [Vec3::ZERO; 2],
        diagonal([0.0, coefficient, 0.0, 0.0, 0.0, 0.0]),
        duration,
        || false,
    )
    .unwrap();
    let result = prepared.relax(states, controls(1e-10), || false).unwrap();
    let after = result.impulse.bodies.map(|receipt| receipt.state_after);
    close(after[0].angular_momentum_body().z, -2.0 * impulse, 2e-12);
    close(after[1].angular_momentum_body().z, 0.0, 2e-12);
    close(after[0].linear_momentum_world().y, -impulse, 2e-12);
    close(after[1].linear_momentum_world().y, 9.0 + impulse, 2e-12);
    close(result.relative_twist_after[1], final_slip, 2e-12);
    for origin in [Vec3::ZERO, Vec3::new(3.0, -4.0, 7.0)] {
        vector(angular(after, origin), angular(states, origin), 2e-11);
    }
    assert!(result.dissipated_energy_j > 0.0);
}

fn full_matrix() -> [[f64; 6]; 6] {
    let coupling = [0.8, -0.4, 0.3, 0.7, -0.2, 0.5];
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            coupling[i] * coupling[j] + if i == j { 0.4 + 0.2 * i as f64 } else { 0.0 }
        })
    })
}
fn full_case() -> ([RigidBodyState; 2], [MassProperties; 2], [Vec3; 2]) {
    let properties = [2.0, 3.0].map(|mass| {
        MassProperties::from_principal_frame(
            mass,
            Vec3::ZERO,
            Vec3::new(2.0, 3.0, 4.0),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_4)
                .unwrap(),
        )
        .unwrap()
    });
    let states = [
        body(
            Vec3::new(-1.0, 2.0, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.4).unwrap(),
            Vec3::new(0.5, -1.0, 2.0),
            Vec3::new(0.8, -0.4, 0.3),
        ),
        body(
            Vec3::new(2.0, -0.5, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.6).unwrap(),
            Vec3::new(-0.7, 0.2, -0.5),
            Vec3::new(0.1, 0.5, -0.6),
        ),
    ];
    (
        states,
        properties,
        [Vec3::new(0.3, -0.2, 0.1), Vec3::new(-0.1, 0.4, 0.2)],
    )
}
fn get(v: &[f64; 13], start: usize) -> Vec3 {
    Vec3::new(v[start], v[start + 1], v[start + 2])
}
fn put(v: &mut [f64; 13], start: usize, value: Vec3) {
    v[start..start + 3].copy_from_slice(&[value.x, value.y, value.z]);
}
fn omega(h: Vec3) -> Vec3 {
    // Independent inverse of [[2.5,-.5,0],[-.5,2.5,0],[0,0,4]].
    Vec3::new(
        (2.5 * h.x + 0.5 * h.y) / 6.0,
        (0.5 * h.x + 2.5 * h.y) / 6.0,
        h.z / 4.0,
    )
}
fn derivative(y: &[f64; 13], poses: [Pose; 2], arms: [Vec3; 2], c: [[f64; 6]; 6]) -> [f64; 13] {
    let r =
        core::array::from_fn::<_, 2, _>(|i| poses[i].orientation().rotate_body_to_world(arms[i]));
    let spin = [
        poses[0]
            .orientation()
            .rotate_body_to_world(omega(get(y, 3))),
        poses[1]
            .orientation()
            .rotate_body_to_world(omega(get(y, 9))),
    ];
    let displacement = poses[1].position_world().sub(poses[0].position_world());
    let slip = get(y, 6)
        .scale(1.0 / 3.0)
        .sub(get(y, 0).scale(1.0 / 2.0))
        .add(spin[1].cross(r[1]))
        .sub(spin[0].cross(displacement.add(r[1])));
    let translation = poses[0].orientation().rotate_world_to_body(slip);
    let rotation = poses[0]
        .orientation()
        .rotate_world_to_body(spin[1].sub(spin[0]));
    let rate = [
        translation.x,
        translation.y,
        translation.z,
        rotation.x,
        rotation.y,
        rotation.z,
    ];
    let effort: [f64; 6] =
        core::array::from_fn(|i| -(0..6).map(|j| c[i][j] * rate[j]).sum::<f64>());
    let force = poses[0]
        .orientation()
        .rotate_body_to_world(Vec3::new(effort[0], effort[1], effort[2]));
    let couple = poses[0]
        .orientation()
        .rotate_body_to_world(Vec3::new(effort[3], effort[4], effort[5]));
    let mut result = [0.0; 13];
    put(&mut result, 0, force.scale(-1.0));
    put(&mut result, 6, force);
    put(
        &mut result,
        3,
        poses[0]
            .orientation()
            .rotate_world_to_body(couple.add(displacement.add(r[1]).cross(force)).scale(-1.0)),
    );
    put(
        &mut result,
        9,
        poses[1]
            .orientation()
            .rotate_world_to_body(couple.add(r[1].cross(force))),
    );
    result[12] = -(0..6).map(|i| effort[i] * rate[i]).sum::<f64>();
    result
}

#[test]
fn all_six_coupled_modes_match_independent_continuous_newton_euler_dynamics() {
    let (states, properties, arms) = full_case();
    let poses = states.map(|s| s.pose());
    let matrix = full_matrix();
    let duration = 0.7;
    let mut reference = [0.0; 13];
    for i in 0..2 {
        put(&mut reference, 6 * i, states[i].linear_momentum_world());
        put(&mut reference, 6 * i + 3, states[i].angular_momentum_body());
    }
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = derivative(&reference, poses, arms, matrix);
        let k2 = derivative(
            &core::array::from_fn(|i| reference[i] + 0.5 * h * k1[i]),
            poses,
            arms,
            matrix,
        );
        let k3 = derivative(
            &core::array::from_fn(|i| reference[i] + 0.5 * h * k2[i]),
            poses,
            arms,
            matrix,
        );
        let k4 = derivative(
            &core::array::from_fn(|i| reference[i] + h * k3[i]),
            poses,
            arms,
            matrix,
        );
        for i in 0..13 {
            reference[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    for scales in [
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
        BushingCoordinateScales {
            translation_m: 0.2,
            rotation_rad: 0.7,
        },
    ] {
        let law = PassiveBushingMatrix::try_new(matrix, scales).unwrap();
        let prepared =
            PreparedViscousConnector::new(poses, properties, arms, law, duration, || false)
                .unwrap();
        let result = prepared.relax(states, controls(1e-9), || false).unwrap();
        for i in 0..2 {
            vector(
                result.impulse.bodies[i].state_after.linear_momentum_world(),
                get(&reference, 6 * i),
                3e-10,
            );
            vector(
                result.impulse.bodies[i].state_after.angular_momentum_body(),
                get(&reference, 6 * i + 3),
                3e-10,
            );
        }
        close(result.dissipated_energy_j, reference[12], 3e-10);
        close(result.kinetic_heat_residual_j, 0.0, 3e-11);
        close(result.connector_heat_residual_j, 0.0, 3e-11);
        vector(
            result.impulse.angular_momentum_change_world_nms,
            Vec3::ZERO,
            3e-11,
        );
    }
}

#[test]
fn common_rigid_motion_is_undamped_even_at_separated_attachments() {
    let properties = scalar_properties();
    let angular_velocity = Vec3::new(0.0, 0.0, 0.5);
    let translation = Vec3::new(3.0, -1.0, 2.0);
    let states = [Vec3::new(-2.0, 1.0, 0.0), Vec3::new(2.0, -1.0, 1.0)].map(|position| {
        let index = usize::from(position.x > 0.0);
        body(
            position,
            UnitQuaternion::IDENTITY,
            translation
                .add(angular_velocity.cross(position))
                .scale(properties[index].mass()),
            Vec3::new(0.0, 0.0, 2.0),
        )
    });
    let prepared = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        [Vec3::new(0.5, -0.25, 0.75), Vec3::new(-0.5, 0.25, 0.5)],
        damping(full_matrix()),
        1000.0,
        || false,
    )
    .unwrap();
    let result = prepared.relax(states, controls(1e-11), || false).unwrap();
    assert_eq!(result.relative_twist_before, [0.0; 6]);
    assert_eq!(result.dissipated_energy_j, 0.0);
    assert_eq!(result.impulse.bodies.map(|r| r.state_after), states);
}

#[test]
fn world_origin_and_common_boost_do_not_change_physical_relaxation() {
    let (states, properties, arms) = full_case();
    let prepared = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        arms,
        damping(full_matrix()),
        0.3,
        || false,
    )
    .unwrap();
    let base = prepared.relax(states, controls(1e-9), || false).unwrap();
    let shift = Vec3::new(2.0_f64.powi(30), -2.0_f64.powi(29), 2.0_f64.powi(28));
    let moved = states.map(|s| {
        body(
            s.pose().position_world().add(shift),
            s.pose().orientation(),
            s.linear_momentum_world(),
            s.angular_momentum_body(),
        )
    });
    let moved_prepared = PreparedViscousConnector::new(
        moved.map(|s| s.pose()),
        properties,
        arms,
        damping(full_matrix()),
        0.3,
        || false,
    )
    .unwrap();
    let mut moved_controls = controls(1e-9);
    moved_controls.impulse.angular_origin_world_m = shift;
    let result = moved_prepared
        .relax(moved, moved_controls, || false)
        .unwrap();
    assert_eq!(result.dissipated_energy_j, base.dissipated_energy_j);
    for i in 0..2 {
        assert_eq!(
            result.impulse.bodies[i].linear_impulse_world,
            base.impulse.bodies[i].linear_impulse_world
        );
        assert_eq!(
            result.impulse.bodies[i].angular_impulse_body,
            base.impulse.bodies[i].angular_impulse_body
        );
    }
    let boost = Vec3::new(1024.0, -512.0, 256.0);
    let boosted = core::array::from_fn(|i| {
        body(
            states[i].pose().position_world(),
            states[i].pose().orientation(),
            states[i]
                .linear_momentum_world()
                .add(boost.scale(properties[i].mass())),
            states[i].angular_momentum_body(),
        )
    });
    let result = prepared.relax(boosted, controls(1e-8), || false).unwrap();
    close(result.dissipated_energy_j, base.dissipated_energy_j, 1e-9);
    for i in 0..2 {
        vector(
            result.impulse.bodies[i].linear_impulse_world,
            base.impulse.bodies[i].linear_impulse_world,
            1e-9,
        );
        vector(
            result.impulse.bodies[i].angular_impulse_body,
            base.impulse.bodies[i].angular_impulse_body,
            1e-9,
        );
    }
}

#[test]
fn zero_law_stale_geometry_cancellation_and_undeliverable_impulses_are_atomic() {
    let (states, properties, arms) = full_case();
    let identity = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        arms,
        diagonal([0.0; 6]),
        1.0,
        || false,
    )
    .unwrap();
    let result = identity.relax(states, controls(1e-12), || false).unwrap();
    assert_eq!(result.impulse.bodies.map(|r| r.state_after), states);
    assert_eq!(result.dissipated_energy_j, 0.0);
    let prepared = PreparedViscousConnector::new(
        states.map(|s| s.pose()),
        properties,
        arms,
        damping(full_matrix()),
        0.3,
        || false,
    )
    .unwrap();
    let mut changed = states;
    changed[0] = body(
        states[0]
            .pose()
            .position_world()
            .add(Vec3::new(0.1, 0.0, 0.0)),
        states[0].pose().orientation(),
        states[0].linear_momentum_world(),
        states[0].angular_momentum_body(),
    );
    assert!(matches!(
        prepared.relax(changed, controls(1e-9), || false),
        Err(ViscousConnectorError::Invalid(_))
    ));
    assert!(
        PreparedViscousConnector::new(
            states.map(|s| s.pose()),
            properties,
            arms,
            damping(full_matrix()),
            -1.0,
            || false
        )
        .is_err()
    );
    for cancel_at in [0, 5, 12] {
        let mut polls = 0;
        let cancelled = prepared.relax(states, controls(1e-9), || {
            let stop = polls == cancel_at;
            polls += 1;
            stop
        });
        assert!(cancelled.is_err());
    }
    let massive_boost = core::array::from_fn(|i| {
        body(
            states[i].pose().position_world(),
            states[i].pose().orientation(),
            Vec3::new(1e30 * properties[i].mass(), 0.0, 0.0),
            states[i].angular_momentum_body(),
        )
    });
    assert!(
        prepared
            .relax(massive_boost, controls(1e-12), || false)
            .is_err()
    );
    // Retrying the original immutable input remains deterministic.
    assert_eq!(
        prepared.relax(states, controls(1e-9), || false).unwrap(),
        prepared.relax(states, controls(1e-9), || false).unwrap(),
    );
}

#[test]
fn cross_coupled_rank_one_free_modes_survive_arbitrarily_stiff_relaxation() {
    let properties = scalar_properties();
    let poses = [Pose::identity(); 2];
    // C = w*w^T with w=(1,1,0,0,0,1). This couples two force axes to
    // a free couple; five independent relative twists are exactly undamped.
    let w = [1.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let law = damping(core::array::from_fn(|i| {
        core::array::from_fn(|j| w[i] * w[j])
    }));
    let prepared =
        PreparedViscousConnector::new(poses, properties, [Vec3::ZERO; 2], law, 1e200, || false)
            .unwrap();
    let free = [
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO),
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0, -3.0, 6.0),
            Vec3::new(4.0, -6.0, 0.0),
        ),
    ];
    let untouched = prepared.relax(free, controls(1e-10), || false).unwrap();
    assert_eq!(untouched.dissipated_energy_j, 0.0);
    assert_eq!(untouched.impulse.bodies.map(|r| r.state_after), free);

    let driven = [
        free[0],
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(9.0, -3.0, 6.0),
            free[1].angular_momentum_body(),
        ),
    ];
    let mobility = 2.0 * (1.0 / 2.0 + 1.0 / 3.0) + 1.0 / 4.0 + 1.0 / 4.0;
    let impulse = -2.0 / mobility;
    let result = prepared.relax(driven, controls(1e-10), || false).unwrap();
    vector(
        result.impulse.bodies[1].linear_impulse_world,
        Vec3::new(impulse, impulse, 0.0),
        2e-12,
    );
    vector(
        result.impulse.bodies[1].angular_impulse_body,
        Vec3::new(0.0, 0.0, impulse),
        2e-12,
    );
    close(result.dissipated_energy_j, 2.0 / mobility, 2e-12);
    close(
        result.relative_twist_after[0]
            + result.relative_twist_after[1]
            + result.relative_twist_after[5],
        0.0,
        2e-12,
    );
    // Free linear and angular axes keep their original stored momenta.
    let after = result.impulse.bodies[1].state_after;
    assert_eq!(after.linear_momentum_world().z, 6.0);
    assert_eq!(after.angular_momentum_body().x, 4.0);
    assert_eq!(after.angular_momentum_body().y, -6.0);
}
