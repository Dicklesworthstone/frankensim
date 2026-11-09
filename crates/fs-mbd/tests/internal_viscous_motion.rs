//! Moving geometry must participate in passive spatial momentum exchange.

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls},
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads, HeldSpatialLoads, InternalWrenchStepError,
    viscous::{
        ViscousMotionError, ViscousMotionLaw, step_viscous_connector,
        step_viscous_connector_with_loads, step_viscous_connector_with_spatial_loads,
    },
};
use fs_mbd::spatial_loads::midpoint::SpatialMidpointConfig;
use fs_mbd::spatial_loads::{PreparedSpatialLoads, WorldPointForce};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn near_vec(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn norm(a: Vec3) -> f64 {
    a.x.hypot(a.y).hypot(a.z)
}
fn controls() -> ViscousConnectorControls {
    ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::new(0.2, -0.3, 0.1),
            linear_tolerance_ns: 2e-10,
            angular_tolerance_nms: 2e-10,
            work_tolerance_j: 2e-10,
        },
        linear_speed_tolerance_m_s: 2e-10,
        angular_speed_tolerance_rad_s: 2e-10,
    }
}
fn matrix(a: [[f64; 6]; 6]) -> PassiveBushingMatrix {
    PassiveBushingMatrix::try_new(
        a,
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
    )
    .unwrap()
}
fn law(a: [[f64; 6]; 6]) -> ViscousMotionLaw {
    ViscousMotionLaw {
        arms_body_m: [Vec3::ZERO; 2],
        damping: matrix(a),
    }
}
fn state(
    p: MassProperties,
    position: Vec3,
    orientation: UnitQuaternion,
    velocity: Vec3,
    omega: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(p.mass()),
        p.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn axial() -> ([RigidBodyState; 2], [MassProperties; 2], ViscousMotionLaw) {
    let properties = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(1.5, 1.5, 1.5)).unwrap(),
    ];
    let states = [
        state(
            properties[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(1.2, 0.0, 0.0),
            Vec3::ZERO,
        ),
        state(
            properties[1],
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(-0.3, 0.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let mut damping = [[0.0; 6]; 6];
    damping[0][0] = 2.3;
    (states, properties, law(damping))
}

#[test]
fn translating_pair_has_exact_recoil_heat_and_second_order_position() {
    let (initial, properties, law) = axial();
    let duration: f64 = 0.7;
    let rate = 2.3 * (1.0 / 2.0 + 1.0 / 3.0);
    let slip0 = -1.5;
    let slip = slip0 * (-rate * duration).exp();
    let center_velocity = (2.0 * 1.2 - 3.0 * 0.3) / 5.0;
    let exact_gap = 2.0 + slip0 * (-(-rate * duration).exp_m1()) / rate;
    let exact_loss = 0.5 * (2.0 * 3.0 / 5.0) * slip0 * slip0 * (-(-2.0 * rate * duration).exp_m1());
    let mut errors = [0.0; 3];
    for (rung, steps) in [1, 2, 4].into_iter().enumerate() {
        let mut current = initial;
        let mut heat = 0.0;
        for _ in 0..steps {
            let result = step_viscous_connector(
                current,
                properties,
                &law,
                duration / steps as f64,
                controls(),
                || false,
            )
            .unwrap();
            heat += result.dissipated_energy_j;
            near(result.work_heat_residual_j, 0.0, 3e-14);
            near_vec(result.linear_momentum_residual_ns, Vec3::ZERO, 2e-14);
            near_vec(result.angular_momentum_residual_nms, Vec3::ZERO, 2e-14);
            current = result.states_after;
        }
        near(
            current[0]
                .center_of_mass_velocity_world(properties[0])
                .unwrap()
                .x,
            center_velocity - 3.0 / 5.0 * slip,
            3e-14,
        );
        near(
            current[1]
                .center_of_mass_velocity_world(properties[1])
                .unwrap()
                .x,
            center_velocity + 2.0 / 5.0 * slip,
            3e-14,
        );
        near(heat, exact_loss, 3e-14);
        near(
            (2.0 * current[0].pose().position_world().x
                + 3.0 * current[1].pose().position_world().x)
                / 5.0,
            1.2 + center_velocity * duration,
            3e-14,
        );
        errors[rung] = (current[1].pose().position_world().x
            - current[0].pose().position_world().x
            - exact_gap)
            .abs();
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.26 * errors[1], "{errors:?}");
}

fn coupled() -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    ViscousMotionLaw,
    [[f64; 6]; 6],
) {
    let properties = [
        MassProperties::from_principal_frame(
            1.7,
            Vec3::ZERO,
            Vec3::new(0.7, 1.1, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.3), 0.4).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            2.4,
            Vec3::ZERO,
            Vec3::new(1.2, 1.7, 2.2),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.1, 0.4), -0.3).unwrap(),
        )
        .unwrap(),
    ];
    let states = [
        state(
            properties[0],
            Vec3::new(-0.4, 0.1, 0.2),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.4, -0.2), 0.3).unwrap(),
            Vec3::new(0.5, -0.2, 0.1),
            Vec3::new(0.4, -0.3, 0.6),
        ),
        state(
            properties[1],
            Vec3::new(0.6, -0.2, -0.1),
            UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.3, 1.0), -0.5).unwrap(),
            Vec3::new(-0.1, 0.3, -0.2),
            Vec3::new(-0.2, 0.5, 0.1),
        ),
    ];
    let mut factor = [[0.0; 6]; 6];
    for i in 0..6 {
        factor[i][i] = 0.6 + 0.1 * i as f64;
    }
    factor[0][3] = 0.25;
    factor[1][5] = -0.2;
    factor[2][4] = 0.15;
    let damping = core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| factor[i][k] * factor[j][k]).sum())
    });
    let mut law = law(damping);
    law.arms_body_m = [Vec3::new(0.2, -0.1, 0.15), Vec3::new(-0.15, 0.3, 0.1)];
    (states, properties, law, damping)
}

// Independent continuous Newton-Euler reference. The response kernel, native
// impulse maps and new splitting are not called from this derivative.
#[derive(Clone, Copy)]
struct Reference {
    x: [Vec3; 2],
    p: [Vec3; 2],
    h: [Vec3; 2],
    q: [[f64; 4]; 2],
    heat: f64,
}
impl Reference {
    fn from_states(states: [RigidBodyState; 2]) -> Self {
        Self {
            x: states.map(|s| s.pose().position_world()),
            p: states.map(|s| s.linear_momentum_world()),
            h: states.map(|s| s.angular_momentum_body()),
            q: states.map(|s| s.pose().orientation().components()),
            heat: 0.0,
        }
    }
    fn shifted(self, r: Self, dt: f64) -> Self {
        Self {
            x: core::array::from_fn(|i| self.x[i].add(r.x[i].scale(dt))),
            p: core::array::from_fn(|i| self.p[i].add(r.p[i].scale(dt))),
            h: core::array::from_fn(|i| self.h[i].add(r.h[i].scale(dt))),
            q: core::array::from_fn(|i| core::array::from_fn(|j| self.q[i][j] + dt * r.q[i][j])),
            heat: self.heat + dt * r.heat,
        }
    }
}
fn derivative(
    s: Reference,
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    c: [[f64; 6]; 6],
) -> Reference {
    let q =
        s.q.map(|q| UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap());
    let omega_body = core::array::from_fn::<_, 2, _>(|i| {
        properties[i].angular_velocity_body_checked(s.h[i]).unwrap()
    });
    let omega = core::array::from_fn::<_, 2, _>(|i| q[i].rotate_body_to_world(omega_body[i]));
    let velocity = core::array::from_fn::<_, 2, _>(|i| s.p[i].scale(1.0 / properties[i].mass()));
    let arms = core::array::from_fn::<_, 2, _>(|i| q[i].rotate_body_to_world(law.arms_body_m[i]));
    let point = core::array::from_fn::<_, 2, _>(|i| s.x[i].add(arms[i]));
    let separation = point[1].sub(point[0]);
    let relative = q[0].rotate_world_to_body(
        velocity[1].add(omega[1].cross(arms[1])).sub(
            velocity[0]
                .add(omega[0].cross(arms[0]))
                .add(omega[0].cross(separation)),
        ),
    );
    let spin = q[0].rotate_world_to_body(omega[1].sub(omega[0]));
    let twist = [relative.x, relative.y, relative.z, spin.x, spin.y, spin.z];
    let load: [f64; 6] = core::array::from_fn(|i| -(0..6).map(|j| c[i][j] * twist[j]).sum::<f64>());
    let force = q[0].rotate_body_to_world(Vec3::new(load[0], load[1], load[2]));
    let couple = q[0].rotate_body_to_world(Vec3::new(load[3], load[4], load[5]));
    let torque = [
        q[0].rotate_world_to_body(
            arms[0]
                .cross(force)
                .scale(-1.0)
                .sub(couple)
                .sub(separation.cross(force)),
        ),
        q[1].rotate_world_to_body(arms[1].cross(force).add(couple)),
    ];
    Reference {
        x: velocity,
        p: [force.scale(-1.0), force],
        h: core::array::from_fn(|i| torque[i].sub(omega_body[i].cross(s.h[i]))),
        q: core::array::from_fn(|i| {
            let [w, x, y, z] = s.q[i];
            let o = omega_body[i];
            [
                -0.5 * (x * o.x + y * o.y + z * o.z),
                0.5 * (w * o.x + y * o.z - z * o.y),
                0.5 * (w * o.y + z * o.x - x * o.z),
                0.5 * (w * o.z + x * o.y - y * o.x),
            ]
        }),
        heat: -(0..6).map(|i| twist[i] * load[i]).sum::<f64>(),
    }
}
fn reference(
    initial: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    c: [[f64; 6]; 6],
    duration: f64,
) -> Reference {
    let mut s = Reference::from_states(initial);
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = derivative(s, properties, law, c);
        let k2 = derivative(s.shifted(k1, 0.5 * dt), properties, law, c);
        let k3 = derivative(s.shifted(k2, 0.5 * dt), properties, law, c);
        let k4 = derivative(s.shifted(k3, dt), properties, law, c);
        s = s
            .shifted(k1, dt / 6.0)
            .shifted(k2, dt / 3.0)
            .shifted(k3, dt / 3.0)
            .shifted(k4, dt / 6.0);
    }
    s
}
fn reference_error(states: [RigidBodyState; 2], heat: f64, reference: Reference) -> f64 {
    let mut error = (heat - reference.heat).abs();
    for i in 0..2 {
        error += norm(states[i].pose().position_world().sub(reference.x[i]));
        error += norm(states[i].linear_momentum_world().sub(reference.p[i]));
        error += norm(states[i].angular_momentum_body().sub(reference.h[i]));
        for (actual, expected) in states[i]
            .pose()
            .orientation()
            .components()
            .into_iter()
            .zip(reference.q[i])
        {
            error += (actual - expected).abs();
        }
    }
    error
}

#[test]
fn full_spatial_moving_geometry_refines_to_independent_newton_euler_flow() {
    let (initial, properties, law, c) = coupled();
    let duration = 0.4;
    let expected = reference(initial, properties, &law, c, duration);
    let mut errors = [0.0; 3];
    for (rung, steps) in [2, 4, 8].into_iter().enumerate() {
        let mut states = initial;
        let mut heat = 0.0;
        for _ in 0..steps {
            let result = step_viscous_connector(
                states,
                properties,
                &law,
                duration / steps as f64,
                controls(),
                || false,
            )
            .unwrap();
            near(result.work_heat_residual_j, 0.0, 2e-11);
            near_vec(result.linear_momentum_residual_ns, Vec3::ZERO, 2e-12);
            near_vec(result.angular_momentum_residual_nms, Vec3::ZERO, 2e-11);
            assert!(result.dissipated_energy_j >= 0.0);
            let endpoint_before = result.half_sources[1]
                .impulse
                .bodies
                .map(|b| b.state_before.pose());
            assert_eq!(endpoint_before, result.states_after.map(|s| s.pose()));
            assert_ne!(endpoint_before, states.map(|s| s.pose()));
            states = result.states_after;
            heat += result.dissipated_energy_j;
        }
        errors[rung] = reference_error(states, heat, expected);
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
}

#[test]
fn zero_law_reproduces_native_loaded_motion_and_gravity_stays_separate_from_heat() {
    let (initial, properties, _) = axial();
    let zero = law([[0.0; 6]; 6]);
    let gravity = Gravity::new(Vec3::new(0.0, -9.81, 0.0)).unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.3, -0.2, 0.1),
            torque_body: Vec3::new(0.1, 0.05, -0.02),
        },
        Wrench {
            force_world: Vec3::new(-0.1, 0.3, 0.0),
            torque_body: Vec3::ZERO,
        },
    ];
    let duration = 0.03;
    let result = step_viscous_connector_with_loads(
        initial,
        properties,
        &zero,
        HeldExternalLoads {
            wrenches: &wrenches,
            gravity,
        },
        duration,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        let native = RigidBodyIntegrator::new(gravity);
        let expected = if wrenches[i].torque_body == Vec3::ZERO {
            native.step(initial[i], properties[i], wrenches[i], duration)
        } else {
            native.step_cayley(initial[i], properties[i], wrenches[i], duration)
        }
        .unwrap();
        assert_eq!(result.states_after[i], expected.state_after);
    }
    assert_eq!(result.dissipated_energy_j, 0.0);
    near(
        result.mechanical_energy_change_j,
        result.external_loads.applied_work_j,
        2e-13,
    );
    near(
        result.gravitational_potential_energy_change_j,
        -result.external_loads.gravity_work_j,
        0.0,
    );

    let (initial, properties, law) = axial();
    let free =
        step_viscous_connector(initial, properties, &law, 0.1, controls(), || false).unwrap();
    let falling = step_viscous_connector_with_loads(
        initial,
        properties,
        &law,
        HeldExternalLoads {
            wrenches: &[Wrench::ZERO; 2],
            gravity,
        },
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    near(free.dissipated_energy_j, falling.dissipated_energy_j, 1e-13);
    for i in 0..2 {
        near_vec(
            falling.states_after[i].linear_momentum_world(),
            free.states_after[i].linear_momentum_world().add(
                gravity
                    .acceleration_world()
                    .scale(0.1 * properties[i].mass()),
            ),
            2e-13,
        );
        near_vec(
            falling.states_after[i].pose().position_world(),
            free.states_after[i]
                .pose()
                .position_world()
                .add(gravity.acceleration_world().scale(0.005)),
            2e-13,
        );
    }
    near(falling.work_heat_residual_j, 0.0, 2e-13);
}

#[test]
fn common_boost_and_world_translation_preserve_spatial_relaxation_and_motion() {
    let (initial, properties, law, _) = coupled();
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let shift = Vec3::new(100.0, -200.0, 50.0);
    let changed = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                initial[i].pose().position_world().add(shift),
                initial[i].pose().orientation(),
            )
            .unwrap(),
            initial[i]
                .linear_momentum_world()
                .add(boost.scale(properties[i].mass())),
            initial[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let mut translated_controls = controls();
    translated_controls.impulse.angular_origin_world_m = translated_controls
        .impulse
        .angular_origin_world_m
        .add(shift);
    let base =
        step_viscous_connector(initial, properties, &law, 0.1, controls(), || false).unwrap();
    let moved = step_viscous_connector(changed, properties, &law, 0.1, translated_controls, || {
        false
    })
    .unwrap();
    near(base.dissipated_energy_j, moved.dissipated_energy_j, 2e-13);
    for i in 0..2 {
        near_vec(
            moved.states_after[i].pose().position_world(),
            base.states_after[i]
                .pose()
                .position_world()
                .add(shift)
                .add(boost.scale(0.1)),
            2e-13,
        );
        near_vec(
            moved.states_after[i].linear_momentum_world(),
            base.states_after[i]
                .linear_momentum_world()
                .add(boost.scale(properties[i].mass())),
            2e-13,
        );
        near_vec(
            moved.states_after[i].angular_momentum_body(),
            base.states_after[i].angular_momentum_body(),
            2e-13,
        );
    }
}

#[test]
fn invalid_and_final_cancellation_refuse_without_consuming_the_original_law() {
    let (initial, properties, law, _) = coupled();
    let saved_law = law.clone();
    let mut polls = 0;
    let valid = step_viscous_connector(initial, properties, &law, 0.1, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut later = 0;
    assert!(matches!(
        step_viscous_connector(initial, properties, &law, 0.1, controls(), || {
            later += 1;
            later == polls
        }),
        Err(ViscousMotionError::Motion(
            InternalWrenchStepError::Cancelled
        ))
    ));
    assert_eq!(later, polls);
    for duration in [0.0, -0.1, f64::NAN, f64::from_bits(1)] {
        assert!(
            step_viscous_connector(initial, properties, &law, duration, controls(), || false)
                .is_err()
        );
    }
    assert!(matches!(
        step_viscous_connector_with_loads(
            initial,
            properties,
            &law,
            HeldExternalLoads {
                wrenches: &[Wrench::ZERO],
                gravity: Gravity::ZERO
            },
            0.1,
            controls(),
            || false
        ),
        Err(ViscousMotionError::Motion(
            InternalWrenchStepError::ExternalLoadCount { .. }
        ))
    ));
    assert_eq!(law, saved_law);
    assert_eq!(
        step_viscous_connector(initial, properties, &law, 0.1, controls(), || false).unwrap(),
        valid
    );
}

#[test]
fn rotating_point_forces_and_world_couples_keep_original_external_work() {
    let (initial, properties, law, _) = coupled();
    let points = [
        WorldPointForce {
            arm_body: Vec3::new(0.2, -0.1, 0.15),
            force_world: Vec3::new(0.1, 0.3, -0.2),
        },
        WorldPointForce {
            arm_body: Vec3::new(-0.3, 0.2, 0.1),
            force_world: Vec3::new(-0.2, 0.1, 0.05),
        },
    ];
    let couples = [Vec3::new(0.03, -0.02, 0.04), Vec3::new(-0.01, 0.03, 0.02)];
    let loads = core::array::from_fn::<_, 2, _>(|i| {
        PreparedSpatialLoads::new(Wrench::ZERO, couples[i], &[points[i]], || false).unwrap()
    });
    let duration = 0.07;
    let result = step_viscous_connector_with_spatial_loads(
        initial,
        properties,
        &law,
        HeldSpatialLoads {
            bodies: &loads,
            gravity: Gravity::ZERO,
            midpoint: SpatialMidpointConfig {
                maximum_iterations: 24,
                maximum_backtracks: 12,
                angular_momentum_tolerance_nms: 1e-13,
                relative_tolerance: 1e-13,
                maximum_cayley_increment: 0.5,
            },
        },
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let mut force_work = 0.0;
    let mut angular_impulse = Vec3::ZERO;
    for i in 0..2 {
        let before = result.half_sources[0].impulse.bodies[i].state_after;
        let after = result.half_sources[1].impulse.bodies[i].state_before;
        let r0 = before
            .pose()
            .orientation()
            .rotate_body_to_world(points[i].arm_body);
        let r1 = after
            .pose()
            .orientation()
            .rotate_body_to_world(points[i].arm_body);
        let x0 = before.pose().position_world().add(r0);
        let x1 = after.pose().position_world().add(r1);
        force_work += points[i].force_world.dot(x1.sub(x0));
        let mean_arm = x0
            .add(x1)
            .scale(0.5)
            .sub(controls().impulse.angular_origin_world_m);
        angular_impulse = angular_impulse.add(
            mean_arm
                .cross(points[i].force_world)
                .add(couples[i])
                .scale(duration),
        );
    }
    near(
        result.external_loads.applied_force_work_j,
        force_work,
        2e-13,
    );
    near_vec(
        result.external_loads.linear_impulse_world_ns,
        points[0]
            .force_world
            .add(points[1].force_world)
            .scale(duration),
        2e-14,
    );
    near_vec(
        result.external_loads.angular_impulse_world_nms,
        angular_impulse,
        2e-13,
    );
    near(
        result.mechanical_energy_change_j + result.dissipated_energy_j,
        result.external_loads.applied_work_j,
        2e-11,
    );
    near(result.work_heat_residual_j, 0.0, 2e-11);
    assert!(result.dissipated_energy_j > 0.0);
}

#[path = "internal_viscous_motion/driven.rs"]
mod driven;
